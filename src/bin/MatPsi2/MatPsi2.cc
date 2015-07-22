
#include "MatPsi2.h"
#include <read_options.cc>

namespace psi {
#ifdef PSIDEBUG
    FILE* outfile = stdout;
#else
    FILE* outfile = fopen("/dev/null", "w");
#endif
    char* psi_file_prefix = "matpsi2";
    std::string outfile_name = "";
    extern int read_options(const std::string &name, Options & options, 
            bool suppress_printing = false);
}

unsigned long int parse_memory_str(const std::string& memory_str) {
    std::string memory_str_ = memory_str;
    boost::algorithm::to_lower(memory_str_);
    boost::cmatch cm;
    boost::regex_search(memory_str_.c_str(), cm, boost::regex("([+]?[0-9]*.?[0-9]+)"));
    double memory = boost::lexical_cast<double>(std::string(cm[1].first, cm[1].second));
    boost::regex_search(memory_str_.c_str(), cm, boost::regex("([a-z])"));
    std::string unit_str(cm[1].first, cm[1].second);
    unsigned long int unit;
    if(boost::iequals(unit_str, "b"))
        unit = 1L;
    else if(boost::iequals(unit_str, "k"))
        unit = 1000L;
    else if(boost::iequals(unit_str, "m"))
        unit = 1000000L;
    else
        unit = 1000000000L;
    if((unsigned long int)(memory*unit) < 100000000L) // less than 100 mb 
        return 100000000L; // if less than 100mb then return 100 mb 
    else
        return (unsigned long int)(memory*unit);
}

// Constructor
MatPsi2::MatPsi2(SharedMatrix cartesian, const std::string& basisname, int charge, int multiplicity, const std::string& path)
{
    // some necessary initializations
    process_environment_.initialize();
    
    // set cores and memory 
    process_environment_.set_n_threads(1);
    process_environment_.set_memory(parse_memory_str("1000mb"));
    worldcomm_ = initialize_communicator(0, NULL, process_environment_);
    process_environment_.set_worldcomm(worldcomm_);
    
    // read in options 
    process_environment_.options.set_read_globals(true);
    read_options("", process_environment_.options, true);
    process_environment_.options.set_read_globals(false);
    process_environment_.set("PSIDATADIR", path);
    process_environment_.options.set_global_int("MAXITER", 100);
    
    Wavefunction::initialize_singletons();
    
    // initialize psio 
    psio_ = boost::shared_ptr<PSIO>(new PSIO);
    process_environment_.set_psio(psio_);
    
    // create molecule object and set its basis set name 
    molecule_ = psi::Molecule::create_molecule_from_cartesian(process_environment_, cartesian, charge, multiplicity);
    molecule_->set_reinterpret_coordentry(false);
    molecule_->set_basis_all_atoms(basisname);
    process_environment_.set_molecule(molecule_);
    
    // create basis object and one & two electron integral factories & rhf 
    create_basis_and_integral_factories();
    
    // create matrix factory object 
    int nbf[] = { basis_->nbf() };
    matfac_ = boost::shared_ptr<MatrixFactory>(new MatrixFactory);
    matfac_->init_with(1, nbf, nbf);
    
    // set default SCF reference according to multiplicity
    if(molecule_->multiplicity() > 1) {
        process_environment_.options.set_global_str("REFERENCE", "UHF");
    } else {
        process_environment_.options.set_global_str("REFERENCE", "RHF");
    }
    
    // set default DFT functional as B3LYP
    process_environment_.options.set_global_str("DFT_FUNCTIONAL", "B3LYP");
}

void MatPsi2::create_basis() {
    // create basis object 
    boost::shared_ptr<PointGroup> c1group(new PointGroup("C1"));
    molecule_->set_point_group(c1group); 
    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    basis_ = BasisSet::construct(process_environment_, parser, molecule_, "BASIS");  
    
    molecule_->set_point_group(c1group); // creating basis set object change molecule's point group, for some reasons 
}

void MatPsi2::create_basis_and_integral_factories() {
    
    create_basis();
    
    // create integral factory object 
    intfac_ = boost::shared_ptr<IntegralFactory>(new IntegralFactory(basis_, basis_, basis_, basis_));
    
    // create two electron integral generator
    eri_ = boost::shared_ptr<TwoBodyAOInt>(intfac_->eri());
}

// destructor 
MatPsi2::~MatPsi2() {
    if(wfn_ != NULL)
        wfn_->extern_finalize();
    if(jk_ != NULL)
        jk_->finalize();
}

void MatPsi2::Settings_SetMaxNumCPUCores(int ncores) {
    // set cores and update worldcomm_ 
    process_environment_.set_n_threads(ncores);
    worldcomm_ = initialize_communicator(0, NULL, process_environment_);
    process_environment_.set_worldcomm(worldcomm_);
}

void MatPsi2::Settings_SetMaxMemory(std::string memory_str) {
    // set memory and update worldcomm_ 
    process_environment_.set_memory(parse_memory_str(memory_str));
    worldcomm_ = initialize_communicator(0, NULL, process_environment_);
    process_environment_.set_worldcomm(worldcomm_);
}

void MatPsi2::Molecule_Fix() {
    molecule_->set_orientation_fixed();
    molecule_->set_com_fixed();
    molecule_->set_reinterpret_coordentry(false);
}

void MatPsi2::Molecule_Free() {
    molecule_->set_orientation_fixed(false);
    molecule_->set_com_fixed(false);
    molecule_->set_reinterpret_coordentry(true);
    molecule_->activate_all_fragments(); // a trick to set lock_frame_ = false;
    molecule_->update_geometry();
}

void MatPsi2::Molecule_SetGeometry(SharedMatrix newGeom) {
    
    // store the old geometry
    Matrix oldgeom = molecule_->geometry();
    molecule_->set_geometry(*(newGeom.get()));
    
    // determine whether the new geometry will cause a problem (typically 2 atoms are at the same point) 
    Matrix distmat = molecule_->distance_matrix();
    for(int i = 0; i < molecule_->natom(); i++) {
        for(int j = 0; j < i - 1; j++) {
            if(distmat.get(i, j) == 0) {
                molecule_->set_geometry(oldgeom);
                throw PSIEXCEPTION("Molecule_SetGeometry: The new geometry has (at least) two atoms at the same spot.");
            }
        }
    }
    
    // update other objects 
    if(jk_ != NULL)
        jk_->finalize();
    psio_->_psio_manager_->psiclean();
    jk_.reset();
    create_basis_and_integral_factories();
}

SharedVector MatPsi2::Molecule_AtomicNumbers() {
    SharedVector zlistvec(new Vector(molecule_->natom()));
    for(int i = 0; i < molecule_->natom(); i++) {
        zlistvec->set(i, (double)molecule_->Z(i));
    }
    return zlistvec;
}

int MatPsi2::Molecule_NumElectrons() {
    int charge = molecule_->molecular_charge();
    int nelectron  = 0;
    for(int i = 0; i < molecule_->natom(); i++)
        nelectron += (int)molecule_->Z(i);
    nelectron -= charge;
    return nelectron;
}

void MatPsi2::Molecule_SetChargeMult(int charge, int mult) {
    SharedVector charge_mult = Molecule_ChargeMult();
    int nelectron = Molecule_NumElectrons() + charge_mult->get(0) - charge;
    if(mult - 1 > nelectron || mult%2 == nelectron%2){
        throw PSIEXCEPTION("Molecule_SetChargeMult: Charge and Multiplicity are not compatible.");
    }
    molecule_->set_molecular_charge(charge);
    molecule_->set_multiplicity(mult);
}

SharedVector MatPsi2::Molecule_ChargeMult() {
    SharedVector charge_mult(new Vector(2));
    charge_mult->set(0, (double)molecule_->molecular_charge());
    charge_mult->set(1, (double)molecule_->multiplicity());
    return charge_mult;
}

void MatPsi2::BasisSet_SetBasisSet(const std::string& basisname) {
    if(jk_ != NULL)
        jk_->finalize();
    psio_->_psio_manager_->psiclean();
    jk_.reset();
    
    molecule_->set_basis_all_atoms(basisname);
    
    // create basis object and one & two electron integral factories & rhf 
    create_basis_and_integral_factories();
    
    // create matrix factory object 
    int nbf[] = { basis_->nbf() };
    matfac_ = boost::shared_ptr<MatrixFactory>(new MatrixFactory);
    matfac_->init_with(1, nbf, nbf);
}

SharedVector MatPsi2::BasisSet_ShellTypes() {
    SharedVector shellTypesVec(new Vector(basis_->nshell()));
    for(int i = 0; i < basis_->nshell(); i++) {
        shellTypesVec->set(i, (double)basis_->shell(i).am());
    }
    return shellTypesVec;
}

SharedVector MatPsi2::BasisSet_ShellNumFunctions() {
    SharedVector shellNfuncsVec(new Vector(basis_->nshell()));
    for(int i = 0; i < basis_->nshell(); i++) {
        shellNfuncsVec->set(i, (double)basis_->shell(i).nfunction());
    }
    return shellNfuncsVec;
}

SharedVector MatPsi2::BasisSet_ShellNumPrimitives() {
    SharedVector shellNprimsVec(new Vector(basis_->nshell()));
    for(int i = 0; i < basis_->nshell(); i++) {
        shellNprimsVec->set(i, (double)basis_->shell(i).nprimitive());
    }
    return shellNprimsVec;
}

SharedVector MatPsi2::BasisSet_ShellToCenter() {
    SharedVector shell2centerVec(new Vector(basis_->nshell()));
    for(int i = 0; i < basis_->nshell(); i++) {
        shell2centerVec->set(i, (double)basis_->shell_to_center(i));
    }
    return shell2centerVec;
}

SharedVector MatPsi2::BasisSet_FuncToCenter() {
    SharedVector func2centerVec(new Vector(basis_->nbf()));
    for(int i = 0; i < basis_->nbf(); i++) {
        func2centerVec->set(i, (double)basis_->function_to_center(i));
    }
    return func2centerVec;
}

SharedVector MatPsi2::BasisSet_FuncToShell() {
    SharedVector func2shellVec(new Vector(basis_->nbf()));
    for(int i = 0; i < basis_->nbf(); i++) {
        func2shellVec->set(i, (double)basis_->function_to_shell(i));
    }
    return func2shellVec;
}

SharedVector MatPsi2::BasisSet_FuncToAngular() {
    SharedVector func2amVec(new Vector(basis_->nbf()));
    for(int i = 0; i < basis_->nbf(); i++) {
        func2amVec->set(i, (double)basis_->shell(basis_->function_to_shell(i)).am());
    }
    return func2amVec;
}

SharedVector MatPsi2::BasisSet_PrimExp() {
    SharedVector primExpsVec(new Vector(basis_->nprimitive()));
    std::vector<double> temp;
    for(int i = 0; i < basis_->nshell(); i++) {
        std::vector<double> currExps = basis_->shell(i).exps();
        temp.insert(temp.end(), currExps.begin(), currExps.end());
    }
    for(int i = 0; i < basis_->nprimitive(); i++) {
        primExpsVec->set(i, temp[i]);
    }
    return primExpsVec;
}

SharedVector MatPsi2::BasisSet_PrimCoeffUnnorm() {
    SharedVector primCoefsVec(new Vector(basis_->nprimitive()));
    std::vector<double> temp;
    for(int i = 0; i < basis_->nshell(); i++) {
        std::vector<double> currCoefs = basis_->shell(i).original_coefs();
        temp.insert(temp.end(), currCoefs.begin(), currCoefs.end());
    }
    for(int i = 0; i < basis_->nprimitive(); i++) {
        primCoefsVec->set(i, temp[i]);
    }
    return primCoefsVec;
}

SharedMatrix MatPsi2::Integrals_Overlap() {
    SharedMatrix sMat(matfac_->create_matrix("Overlap"));
    boost::shared_ptr<OneBodyAOInt> sOBI(intfac_->ao_overlap());
    sOBI->compute(sMat);
    sMat->hermitivitize();
    return sMat;
}

SharedMatrix MatPsi2::Integrals_Kinetic() {
    SharedMatrix tMat(matfac_->create_matrix("Kinetic"));
    boost::shared_ptr<OneBodyAOInt> tOBI(intfac_->ao_kinetic());
    tOBI->compute(tMat);
    tMat->hermitivitize();
    return tMat;
}

SharedMatrix MatPsi2::Integrals_Potential() {
    SharedMatrix vMat(matfac_->create_matrix("Potential"));
    boost::shared_ptr<OneBodyAOInt> vOBI(intfac_->ao_potential());
    vOBI->compute(vMat);
    vMat->hermitivitize();
    return vMat;
}

std::vector<SharedMatrix> MatPsi2::Integrals_Dipole() {
    std::vector<SharedMatrix> ao_dipole;
    SharedMatrix dipole_x(matfac_->create_matrix("Dipole x"));
    SharedMatrix dipole_y(matfac_->create_matrix("Dipole y"));
    SharedMatrix dipole_z(matfac_->create_matrix("Dipole z"));
    ao_dipole.push_back(dipole_x);
    ao_dipole.push_back(dipole_y);
    ao_dipole.push_back(dipole_z);
    boost::shared_ptr<OneBodyAOInt> dipoleOBI(intfac_->ao_dipole());
    dipoleOBI->compute(ao_dipole);
    ao_dipole[0]->hermitivitize();
    ao_dipole[1]->hermitivitize();
    ao_dipole[2]->hermitivitize();
    return ao_dipole;
}

std::vector<SharedMatrix> MatPsi2::Integrals_PotentialEachCore() {
    int natom = molecule_->natom();
    std::vector<SharedMatrix> viMatVec;
    boost::shared_ptr<OneBodyAOInt> viOBI(intfac_->ao_potential());
    boost::shared_ptr<PotentialInt> viPtI = boost::static_pointer_cast<PotentialInt>(viOBI);
    SharedMatrix Zxyz = viPtI->charge_field();
    SharedMatrix Zxyz_rowi(new Matrix(1, 4));
    for( int i = 0; i < natom; i++) {
        SharedVector Zxyz_rowi_vec = Zxyz->get_row(0, i);
        Zxyz_rowi->set_row(0, 0, Zxyz_rowi_vec);
        viPtI->set_charge_field(Zxyz_rowi);
        viMatVec.push_back(matfac_->create_shared_matrix("PotentialEachCore"));
        viOBI->compute(viMatVec[i]);
        viMatVec[i]->hermitivitize();
    }
    return viMatVec;
}

SharedMatrix MatPsi2::Integrals_PotentialPtQ(SharedMatrix Zxyz_list) {
    boost::shared_ptr<OneBodyAOInt> viOBI(intfac_->ao_potential());
    boost::shared_ptr<PotentialInt> viPtI = boost::static_pointer_cast<PotentialInt>(viOBI);
    viPtI->set_charge_field(Zxyz_list);
    SharedMatrix vZxyzListMat(matfac_->create_matrix("PotentialPointCharges"));
    viOBI->compute(vZxyzListMat);
    vZxyzListMat->hermitivitize();
    return vZxyzListMat;
}

double MatPsi2::Integrals_ijkl(int i, int j, int k, int l) {
    int ish = basis_->function_to_shell(i);
    int jsh = basis_->function_to_shell(j);
    int ksh = basis_->function_to_shell(k);
    int lsh = basis_->function_to_shell(l);
    int ii = i - basis_->shell_to_basis_function(ish);
    int jj = j - basis_->shell_to_basis_function(jsh);
    int kk = k - basis_->shell_to_basis_function(ksh);
    int ll = l - basis_->shell_to_basis_function(lsh);
    int ni = basis_->shell(ish).nfunction();
    int nj = basis_->shell(jsh).nfunction();
    int nk = basis_->shell(ksh).nfunction();
    int nl = basis_->shell(lsh).nfunction();
    eri_->compute_shell(ish, jsh, ksh, lsh);
    const double *buffer = eri_->buffer();
    return buffer[ll+nl*(kk+nk*(jj+nj*ii))];
}

inline int ij2I(int i, int j) {
    if(i < j) {
        int tmp = i;
        i = j;
        j = tmp;
    }
    return i * ( i + 1 ) / 2 + j;
}

int MatPsi2::Integrals_NumUniqueTEIs() {
    int nbf = basis_->nbf();
    return ( nbf * ( nbf + 1 ) * ( nbf * nbf + nbf + 2 ) ) / 8;
}

void MatPsi2::Integrals_AllUniqueTEIs(double* matpt) {
    AOShellCombinationsIterator shellIter = intfac_->shells_iterator();
    int nuniq = Integrals_NumUniqueTEIs();
    const double *buffer = eri_->buffer();
    for (shellIter.first(); shellIter.is_done() == false; shellIter.next()) {
        // Compute quartet
        eri_->compute_shell(shellIter);
        // From the quartet get all the integrals
        AOIntegralsIterator intIter = shellIter.integrals_iterator();
        for (intIter.first(); intIter.is_done() == false; intIter.next()) {
            matpt[ ij2I( ij2I(intIter.i(), intIter.j()), ij2I(intIter.k(), intIter.l()) ) ] = buffer[intIter.index()];
        }
    }
}

void MatPsi2::Integrals_AllTEIs(double* matpt) {
    AOShellCombinationsIterator shellIter = intfac_->shells_iterator();
    int nbf = basis_->nbf();
    const double *buffer = eri_->buffer();
    for (shellIter.first(); shellIter.is_done() == false; shellIter.next()) {
        // Compute quartet
        eri_->compute_shell(shellIter);
        // From the quartet get all the integrals
        AOIntegralsIterator intIter = shellIter.integrals_iterator();
        for (intIter.first(); intIter.is_done() == false; intIter.next()) {
            int i = intIter.i();
            int j = intIter.j();
            int k = intIter.k();
            int l = intIter.l();
            matpt[ l+nbf*(k+nbf*(j+nbf*i)) ] = buffer[intIter.index()];
            matpt[ l+nbf*(k+nbf*(i+nbf*j)) ] = buffer[intIter.index()];
            matpt[ k+nbf*(l+nbf*(j+nbf*i)) ] = buffer[intIter.index()];
            matpt[ k+nbf*(l+nbf*(i+nbf*j)) ] = buffer[intIter.index()];
            matpt[ j+nbf*(i+nbf*(l+nbf*k)) ] = buffer[intIter.index()];
            matpt[ j+nbf*(i+nbf*(k+nbf*l)) ] = buffer[intIter.index()];
            matpt[ i+nbf*(j+nbf*(l+nbf*k)) ] = buffer[intIter.index()];
            matpt[ i+nbf*(j+nbf*(k+nbf*l)) ] = buffer[intIter.index()];
        }
    }
}

void MatPsi2::Integrals_IndicesForK(double* indices1, double* indices2) {
    AOShellCombinationsIterator shellIter = intfac_->shells_iterator();
    for (shellIter.first(); shellIter.is_done() == false; shellIter.next()) {
        // Compute quartet
        //~ eri_->compute_shell(shellIter);
        // From the quartet get all the integrals
        AOIntegralsIterator intIter = shellIter.integrals_iterator();
        for (intIter.first(); intIter.is_done() == false; intIter.next()) {
            int i = intIter.i();
            int j = intIter.j();
            int k = intIter.k();
            int l = intIter.l();
            indices1[ij2I( ij2I(i, j), ij2I(k, l) )] = ij2I( ij2I(i, l), ij2I(k, j) );
            indices2[ij2I( ij2I(i, j), ij2I(k, l) )] = ij2I( ij2I(i, k), ij2I(j, l) );
        }
    }
}

void MatPsi2::JK_Initialize(std::string jktype, std::string auxBasisName) {
    std::transform(jktype.begin(), jktype.end(), jktype.begin(), ::toupper);
    if(jk_ != NULL)
        jk_->finalize();
    if(wfn_ == NULL) {
        std::string scfType = process_environment_.options.get_str("REFERENCE");
        if(scfType == "RHF" || scfType == "RKS") {
            if(molecule_->multiplicity() > 1)
                throw PSIEXCEPTION("JK_Initialize: RHF or RKS can handle singlets only.");
            wfn_ = boost::shared_ptr<scf::RHF>(new scf::RHF(process_environment_, basis_));
        } else if(scfType == "UHF" || scfType == "UKS") {
            wfn_ = boost::shared_ptr<scf::UHF>(new scf::UHF(process_environment_, basis_));
        } else {
            throw PSIEXCEPTION("JK_Initialize: Reference SCF type not recognized.");
        }
        process_environment_.set_wavefunction(wfn_);
        wfn_->extern_finalize();
        wfn_.reset();
    }
    if(jktype == "PKJK") {
        jk_ = boost::shared_ptr<JK>(new PKJK(process_environment_, basis_, psio_));
    } else if(jktype == "DFJK") {
        boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
        molecule_->set_basis_all_atoms(auxBasisName, "DF_BASIS_SCF");
        boost::shared_ptr<BasisSet> auxiliary = BasisSet::construct(process_environment_, parser, molecule_, "DF_BASIS_SCF");
        jk_ = boost::shared_ptr<JK>(new DFJK(process_environment_, basis_, auxiliary, psio_));
        molecule_->set_basis_all_atoms(basis_->name());
    } else if(jktype == "ICJK") {
        jk_ = boost::shared_ptr<JK>(new ICJK(process_environment_, basis_));
    } else if(jktype == "DIRECTJK") {
        jk_ = boost::shared_ptr<JK>(new DirectJK(process_environment_, basis_));
    } else {
        throw PSIEXCEPTION("JK_Initialize: JK type not recognized.");
    }
    jk_->set_memory(process_environment_.get_memory());
    jk_->set_cutoff(0.0);
    jk_->initialize();
}

const std::string& MatPsi2::JK_Type() {
    if(jk_ == NULL)
        throw PSIEXCEPTION("JK_Type: JK object has not been initialized.");
    return jk_->JKtype();
}

SharedMatrix DensToEigVectors(SharedMatrix density) {
    if(density == NULL) {
        return density;
    }
    int dim = density->ncol();
    SharedMatrix eigVectors(new Matrix(dim, dim));
    boost::shared_ptr<Vector> eigValues(new Vector(dim));
    density->diagonalize(eigVectors, eigValues);
    for(int i = 0; i < dim; i++) {
        eigValues->set(i, sqrt(max(0.0, eigValues->get(i))));
    }
    for(int i = 0; i < dim; i++) {
        for(int j = 0; j < dim; j++) {
            eigVectors->set(i, j, eigVectors->get(i, j) * eigValues->get(j));
        }
    }
    return eigVectors;
}

std::vector<SharedMatrix> MatPsi2::JK_DensToJ(SharedMatrix densAlpha, SharedMatrix densBeta) {
	return JK_OccOrbToJ(DensToEigVectors(densAlpha), DensToEigVectors(densBeta)); 
}

std::vector<SharedMatrix> MatPsi2::JK_DensToK(SharedMatrix densAlpha, SharedMatrix densBeta) {
	return JK_OccOrbToK(DensToEigVectors(densAlpha), DensToEigVectors(densBeta)); 
}

std::vector<SharedMatrix> MatPsi2::JK_OccOrbToJ(SharedMatrix occOrbAlpha, SharedMatrix occOrbBeta) {
    if(jk_ == NULL)
        JK_Initialize("PKJK");
    jk_->set_do_K(false);
    JK_CalcAllFromOccOrb(occOrbAlpha, occOrbBeta);
    jk_->set_do_K(true);
    return jk_->J();
}

std::vector<SharedMatrix> MatPsi2::JK_OccOrbToK(SharedMatrix occOrbAlpha, SharedMatrix occOrbBeta) {
    if(jk_ == NULL)
        JK_Initialize("PKJK");
    jk_->set_do_J(false);
    JK_CalcAllFromOccOrb(occOrbAlpha, occOrbBeta);
    jk_->set_do_J(true);
    return jk_->K();
}

void MatPsi2::JK_CalcAllFromDens(SharedMatrix densAlpha, SharedMatrix densBeta) {
    JK_CalcAllFromOccOrb(DensToEigVectors(densAlpha), DensToEigVectors(densBeta));
}

void MatPsi2::JK_CalcAllFromOccOrb(SharedMatrix occOrbAlpha, SharedMatrix occOrbBeta) {
    if(jk_ == NULL)
        JK_Initialize("PKJK");
    jk_->C_left().clear();
    jk_->C_left().push_back(occOrbAlpha);
    if(occOrbBeta != NULL)
        jk_->C_left().push_back(occOrbBeta);
    jk_->compute();
}

std::vector<SharedMatrix> MatPsi2::JK_RetrieveJ() {
	if(jk_ == NULL)
		throw PSIEXCEPTION("JK_RetriveJ: J/K calculation has not been done.");
    return jk_->J();
}

std::vector<SharedMatrix> MatPsi2::JK_RetrieveK() {
	if(jk_ == NULL)
		throw PSIEXCEPTION("JK_RetriveK: J/K calculation has not been done.");
    return jk_->K();
}

void MatPsi2::jk_DFException(std::string functionName) {
    if(jk_ == NULL) {
        JK_Initialize("DFJK");
    }
    if((!boost::iequals(jk_->JKtype(), "DFJK")) || dynamic_cast<DFJK*>(jk_.get()) == NULL) {
        throw PSIEXCEPTION(functionName + ": Can only be used with Density-fitting JK.");
    }
}

SharedMatrix MatPsi2::JK_DFTensor_AuxPriPairs() {
    jk_DFException("JK_DFTensor_AuxPriPairs");
    return boost::static_pointer_cast<DFJK>(jk_)->GetQmn();
}

std::vector<SharedMatrix> MatPsi2::JK_DFTensor_AuxPriPri() {
    jk_DFException("JK_DFTensor_AuxPriPri");
    SharedMatrix QmnUnique = boost::static_pointer_cast<DFJK>(jk_)->GetQmn();
    std::vector<SharedMatrix> QmnFull;
    for(int Q = 0; Q < QmnUnique->nrow(); Q++) {
        QmnFull.push_back(SharedMatrix(new Matrix(basis_->nbf(), basis_->nbf())));
        QmnFull[Q]->set(QmnUnique->const_pointer()[Q]);
    }
    return QmnFull;
}

SharedMatrix MatPsi2::JK_DFMetric_InvJHalf() {
    jk_DFException("JK_DFMetric_InvJHalf");
    return boost::static_pointer_cast<DFJK>(jk_)->GetInvJHalf();
}

void MatPsi2::DFT_Initialize(std::string functionalName) {
    std::transform(functionalName.begin(), functionalName.end(), functionalName.begin(), ::toupper);
    process_environment_.options.set_global_str("DFT_FUNCTIONAL", functionalName);
    std::string scfType = process_environment_.options.get_str("REFERENCE");
    std::string dftType;
    if(scfType == "RHF" || scfType == "RKS") {
        if(molecule_->multiplicity() > 1)
            throw PSIEXCEPTION("DFT_Initialize: RHF or RKS can handle singlets only.");
        dftType = "RV";
    } else if(scfType == "UHF" || scfType == "UKS") {
        dftType = "UV";
    } else {
        throw PSIEXCEPTION("DFT_Initialize: Reference SCF type not recognized.");
    }
    dftPotential_ = VBase::build_V(process_environment_, process_environment_.options, dftType);
    dftPotential_->initialize();
}

std::vector<SharedMatrix> MatPsi2::DFT_DensToV(SharedMatrix densAlpha, SharedMatrix densBeta) {
    return DFT_OccOrbToV(DensToEigVectors(densAlpha), DensToEigVectors(densBeta));
}

std::vector<SharedMatrix> MatPsi2::DFT_OccOrbToV(SharedMatrix occOrbAlpha, SharedMatrix occOrbBeta) {
    if(dftPotential_ == NULL)
        DFT_Initialize(process_environment_.options.get_str("DFT_FUNCTIONAL"));
    if(dynamic_cast<UV*>(dftPotential_.get()) != NULL && occOrbBeta == NULL)
        throw PSIEXCEPTION("DFT_OccOrbToV: Unrestricted functional requires both alpha and beta orbital.");
    std::vector<SharedMatrix> & C = dftPotential_->C();
    C.clear();
    C.push_back(occOrbAlpha);
    if(occOrbBeta != NULL)
        C.push_back(occOrbBeta);
    
    // Run the potential object
    dftPotential_->compute();
    return dftPotential_->V();
}

double MatPsi2::DFT_EnergyXC() {
    std::map<std::string, double>& quad = dftPotential_->quadrature_values();  
    return quad["FUNCTIONAL"];
}

void MatPsi2::SCF_SetSCFType(std::string scfType) {
    std::transform(scfType.begin(), scfType.end(), scfType.begin(), ::toupper);
    process_environment_.options.set_global_str("REFERENCE", scfType);
    process_environment_.options.set_global_str("GUESS", "CORE");
}

void MatPsi2::SCF_SetGuessOrb(SharedMatrix guessOrb) {
    process_environment_.options.set_global_str("GUESS", "ORBITAL");
    guessOrbital_ = guessOrb;
}

void MatPsi2::create_wfn() {
    if(jk_ == NULL)
        JK_Initialize("PKJK");
    std::string scfType = process_environment_.options.get_str("REFERENCE");
    if(scfType == "RHF") {
        if(molecule_->multiplicity() > 1)
            throw PSIEXCEPTION("create_wfn: RHF can handle singlets only.");
        wfn_ = boost::shared_ptr<scf::RHF>(new scf::RHF(process_environment_, jk_));
    } else if(scfType == "UHF") {
        wfn_ = boost::shared_ptr<scf::UHF>(new scf::UHF(process_environment_, jk_));
    } else if(scfType == "RKS") {
        if(molecule_->multiplicity() > 1)
            throw PSIEXCEPTION("create_wfn: RKS can handle singlets only.");
        wfn_ = boost::shared_ptr<scf::RKS>(new scf::RKS(process_environment_, jk_));
    } else if(scfType == "UKS") {
        wfn_ = boost::shared_ptr<scf::UKS>(new scf::UKS(process_environment_, jk_));
    } else {
        throw PSIEXCEPTION("create_wfn: SCF type not recognized.");
    }
}

double MatPsi2::SCF_RunSCF() {
    create_wfn();
    if(process_environment_.options.get_str("GUESS") == "ORBITAL") {
        wfn_->SetGuessOrbital(guessOrbital_);
    }
    process_environment_.set_wavefunction(wfn_);
    return wfn_->compute_energy();
}

void MatPsi2::SCF_EnableMOM(int mom_start) {
    process_environment_.options.set_global_int("MOM_START", mom_start);
}

void MatPsi2::SCF_EnableDamping(double damping_percentage) {
    process_environment_.options.set_global_double("DAMPING_PERCENTAGE", damping_percentage);
}

void MatPsi2::SCF_DisableDIIS() {
    process_environment_.options.set_global_bool("DIIS", false);
    process_environment_.options.set_global_int("MAXITER", 500);
}

void MatPsi2::SCF_EnableDIIS() {
    process_environment_.options.set_global_int("DIIS", true);
    process_environment_.options.set_global_int("MAXITER", 100);
}

void MatPsi2::SCF_GuessSAD() {
    process_environment_.options.set_global_str("GUESS", "SAD");
}

void MatPsi2::SCF_GuessCore() {
    process_environment_.options.set_global_str("GUESS", "CORE");
}

double MatPsi2::SCF_TotalEnergy() { 
    if(wfn_ == NULL) {
        SCF_RunSCF();
    }
    return wfn_->EHF(); 
}

SharedMatrix MatPsi2::SCF_OrbitalAlpha() { 
    if(wfn_ == NULL) {
        SCF_RunSCF();
    }
    return wfn_->Ca(); 
}

SharedMatrix MatPsi2::SCF_OrbitalBeta() { 
    if(wfn_ == NULL) {
        SCF_RunSCF();
    }
    return wfn_->Cb(); 
}

SharedVector MatPsi2::SCF_OrbEigValAlpha() { 
    if(wfn_ == NULL) {
        SCF_RunSCF();
    }
    return wfn_->epsilon_a(); 
}

SharedVector MatPsi2::SCF_OrbEigValBeta() { 
    if(wfn_ == NULL) {
        SCF_RunSCF();
    }
    return wfn_->epsilon_b(); 
}

SharedMatrix MatPsi2::SCF_DensityAlpha() { 
    if(wfn_ == NULL) {
        SCF_RunSCF();
    }
    return wfn_->Da(); 
}

SharedMatrix MatPsi2::SCF_DensityBeta() { 
    if(wfn_ == NULL) {
        SCF_RunSCF();
    }
    return wfn_->Db(); 
}

SharedMatrix MatPsi2::SCF_CoreHamiltonian() {
    if(wfn_ == NULL) {
        SCF_RunSCF();
    }
    return wfn_->H(); 
}

SharedMatrix MatPsi2::SCF_FockAlpha() { 
    if(wfn_ == NULL) {
        SCF_RunSCF();
    }
    return wfn_->Fa(); 
}

SharedMatrix MatPsi2::SCF_FockBeta() { 
    if(wfn_ == NULL) {
        SCF_RunSCF();
    }
    return wfn_->Fb(); 
}

SharedMatrix MatPsi2::SCF_Gradient() {
    if(wfn_ == NULL) {
        SCF_RunSCF();
    }
    scfgrad::SCFGrad scfgrad_ = scfgrad::SCFGrad(process_environment_);
    return scfgrad_.compute_gradient();
}

SharedMatrix MatPsi2::SCF_GuessDensity() {
    create_wfn();
    process_environment_.set_wavefunction(wfn_);
    return wfn_->GuessDensity();
}

SharedMatrix MatPsi2::SCF_RHF_J() { 
    if(wfn_ == NULL) {
        SCF_RunSCF();
    }
    if(dynamic_cast<scf::RHF*>(wfn_.get()) == NULL)
        throw PSIEXCEPTION("SCF_RHF_J: This function works only for RHF.");
    return boost::static_pointer_cast<scf::RHF>(wfn_)->J(); 
}

SharedMatrix MatPsi2::SCF_RHF_K() { 
    if(wfn_ == NULL) {
        SCF_RunSCF();
    }
    if(dynamic_cast<scf::RHF*>(wfn_.get()) == NULL)
        throw PSIEXCEPTION("SCF_RHF_K: This function works only for RHF.");
    return boost::static_pointer_cast<scf::RHF>(wfn_)->K(); 
}


