
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
MatPsi2::MatPsi2(const std::string& molstring, const std::string& basisname, int charge, int multiplicity, const std::string& path)
    : molstring_(molstring), basisname_(basisname)
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
    boost::filesystem::path uniqname = boost::filesystem::unique_path();
    std::string matpsi_id = uniqname.string();
    boost::filesystem::path tempdir = boost::filesystem::temp_directory_path();
    std::string matpsi_tempdir_str = tempdir.string();
    matpsi_tempdir_str += "/matpsi2.temp.";
    matpsi_tempdir_str += matpsi_id;
    psio_ = boost::shared_ptr<PSIO>(new PSIO);
    psio_->set_pid(matpsi_id);
    for (int i=1; i<=PSIO_MAXVOL; ++i) {
        char kwd[20];
        sprintf(kwd, "VOLUME%u", i);
        psio_->filecfg_kwd("DEFAULT", kwd, PSIF_CHKPT, matpsi_tempdir_str.c_str());
        psio_->filecfg_kwd("DEFAULT", kwd, -1, matpsi_tempdir_str.c_str());
    }
    psio_->_psio_manager_->set_default_path(matpsi_tempdir_str);
    process_environment_.set_psio(psio_);
    
    // create molecule object and set its basis set name 
    molecule_ = psi::Molecule::create_molecule_from_string(process_environment_, molstring_);
    molecule_->set_molecular_charge(charge);
    molecule_->set_multiplicity(multiplicity);
    molecule_->set_com_fixed();
    molecule_->set_orientation_fixed();
    molecule_->update_geometry();
    molecule_->set_reinterpret_coordentry(false);
    molecule_->set_basis_all_atoms(basisname_);
    process_environment_.set_molecule(molecule_);
    
    // create basis object and one & two electron integral factories & rhf 
    create_basis_and_integral_factories();
    
    // create matrix factory object 
    int nbf[] = { basis_->nbf() };
    matfac_ = boost::shared_ptr<MatrixFactory>(new MatrixFactory);
    matfac_->init_with(1, nbf, nbf);
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
    if(rhf_ != NULL)
        rhf_->extern_finalize();
    if(jk_ != NULL)
        jk_->finalize();
    psio_->_psio_manager_->psiclean();
    boost::filesystem::remove_all(psio_->_psio_manager_->get_file_path(0));
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

void MatPsi2::BasisSet_SetBasisSet(const std::string& basisname) {
    if(jk_ != NULL)
        jk_->finalize();
    psio_->_psio_manager_->psiclean();
    jk_.reset();
    
    basisname_ = basisname;
    molecule_->set_basis_all_atoms(basisname_);
    
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

SharedVector MatPsi2::BasisSet_FunctionToCenter() {
    SharedVector func2centerVec(new Vector(basis_->nbf()));
    for(int i = 0; i < basis_->nbf(); i++) {
        func2centerVec->set(i, (double)basis_->function_to_center(i));
    }
    return func2centerVec;
}

SharedVector MatPsi2::BasisSet_FunctionToAngularMomentum() {
    SharedVector func2amVec(new Vector(basis_->nbf()));
    for(int i = 0; i < basis_->nbf(); i++) {
        func2amVec->set(i, (double)basis_->shell(basis_->function_to_shell(i)).am());
    }
    return func2amVec;
}

SharedVector MatPsi2::BasisSet_PrimitiveExponents() {
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

SharedVector MatPsi2::BasisSet_PrimitiveCoefficients() { // should have been normalized 
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

SharedMatrix MatPsi2::Integrals_PotentialPointCharges(SharedMatrix Zxyz_list) {
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

void MatPsi2::Integrals_IndicesForExchange(double* indices1, double* indices2) {
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

void MatPsi2::JK_Initialize(std::string jktype, std::string auxiliaryBasisSetName) {
    if(jk_ != NULL)
        jk_->finalize();
    if(process_environment_.wavefunction() == NULL)
        RHF_Reset();
    if(boost::iequals(jktype, "PKJK")) {
        jk_ = boost::shared_ptr<JK>(new PKJK(process_environment_, basis_, psio_));
    } else if(boost::iequals(jktype, "DFJK")) {
        boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
        molecule_->set_basis_all_atoms(auxiliaryBasisSetName, "DF_BASIS_SCF");
        boost::shared_ptr<BasisSet> auxiliary = BasisSet::construct(process_environment_, parser, molecule_, "DF_BASIS_SCF");
        jk_ = boost::shared_ptr<JK>(new DFJK(process_environment_, basis_, auxiliary, psio_));
        molecule_->set_basis_all_atoms(basisname_);
    } else if(boost::iequals(jktype, "ICJK")) {
        jk_ = boost::shared_ptr<JK>(new ICJK(process_environment_, basis_));
    } else if(boost::iequals(jktype, "DirectJK")) {
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

SharedMatrix MatPsi2::JK_DensityToJ(SharedMatrix density) {
    if(jk_ == NULL) {
        JK_Initialize("PKJK");
    }
    jk_->set_do_K(false);
    jk_->C_left().clear();
    jk_->D().clear();
    jk_->D().push_back(density);
    
    jk_->compute_from_D();
    SharedMatrix Jnew = jk_->J()[0];
    Jnew->hermitivitize();
    jk_->set_do_K(true);
    return Jnew;
}

SharedMatrix MatPsi2::JK_DensityToK(SharedMatrix density) {
    if(jk_ == NULL) {
        JK_Initialize("PKJK");
    }
    if(boost::iequals(jk_->JKtype(), "DFJK")) {
        throw PSIEXCEPTION("JK_DensityToK: DensityToK cannot be used with DFJK.");
    }
    jk_->set_do_J(false);
    jk_->C_left().clear();
    jk_->D().clear();
    jk_->D().push_back(density);
    
    jk_->compute_from_D();
    SharedMatrix Knew = jk_->K()[0];
    Knew->hermitivitize();
    jk_->set_do_J(true);
    return Knew;
}

SharedMatrix MatPsi2::JK_OrbitalToJ(SharedMatrix orbital) {
    SharedMatrix occupiedOrbital(new Matrix(orbital->nrow(), Molecule_NumElectrons()/2));
    for(int i = 0; i < Molecule_NumElectrons()/2; i++) {
        occupiedOrbital->set_column(0, i, orbital->get_column(0, i));
    }
    return JK_OccupiedOrbitalToJ(occupiedOrbital);
}

SharedMatrix MatPsi2::JK_OrbitalToK(SharedMatrix orbital) {
    SharedMatrix occupiedOrbital(new Matrix(orbital->nrow(), Molecule_NumElectrons()/2));
    for(int i = 0; i < Molecule_NumElectrons()/2; i++) {
        occupiedOrbital->set_column(0, i, orbital->get_column(0, i));
    }
    return JK_OccupiedOrbitalToK(occupiedOrbital);
}

SharedMatrix MatPsi2::JK_OccupiedOrbitalToJ(SharedMatrix occupiedOrbital) {
    if(jk_ == NULL) {
        JK_Initialize("PKJK");
    }
    jk_->set_do_K(false);
    jk_->C_left().clear();
    jk_->C_left().push_back(occupiedOrbital);
    
    jk_->compute();
    SharedMatrix Jnew = jk_->J()[0];
    Jnew->hermitivitize();
    jk_->set_do_K(true);
    return Jnew;
}

SharedMatrix MatPsi2::JK_OccupiedOrbitalToK(SharedMatrix occupiedOrbital) {
    if(jk_ == NULL) {
        JK_Initialize("PKJK");
    }
    jk_->set_do_J(false);
    jk_->C_left().clear();
    jk_->C_left().push_back(occupiedOrbital);
    
    jk_->compute();
    SharedMatrix Knew = jk_->K()[0];
    Knew->hermitivitize();
    jk_->set_do_J(true);
    return Knew;
}

void MatPsi2::DFJKException(std::string functionName) {
    if(jk_ == NULL) {
        JK_Initialize("DFJK");
    }
    if(!boost::iequals(jk_->JKtype(), "DFJK")) {
        throw PSIEXCEPTION(functionName + ": Can only be used with DFJK.");
    }
    if(!boost::static_pointer_cast<DFJK>(jk_)->IsCore()) {
        throw PSIEXCEPTION(functionName + ": Can only be used with DFJK.");
    }
}

SharedMatrix MatPsi2::DFJK_mnQMatrixUnique() {
    DFJKException("DFJK_mnQMatrixUnique");
    return boost::static_pointer_cast<DFJK>(jk_)->GetQmn()->transpose();
}

std::vector<SharedMatrix> MatPsi2::DFJK_mnQTensorFull() {
    DFJKException("DFJK_mnQTensorFull");
    SharedMatrix QmnUnique = boost::static_pointer_cast<DFJK>(jk_)->GetQmn();
    std::vector<SharedMatrix> mnQFull;
    for(int Q = 0; Q < QmnUnique->nrow(); Q++) {
        mnQFull.push_back(SharedMatrix(new Matrix(basis_->nbf(), basis_->nbf())));
        mnQFull[Q]->set(QmnUnique->const_pointer()[Q]);
    }
    return mnQFull;
}

SharedMatrix MatPsi2::DFJK_mnAMatrixUnique() {
    DFJKException("DFJK_mnAMatrixUnique");
    return boost::static_pointer_cast<DFJK>(jk_)->GetAmn()->transpose();
}

std::vector<SharedMatrix> MatPsi2::DFJK_mnATensorFull() {
    DFJKException("DFJK_mnATensorFull");
    SharedMatrix AmnUnique = boost::static_pointer_cast<DFJK>(jk_)->GetAmn();
    std::vector<SharedMatrix> mnAFull;
    for(int aux = 0; aux < AmnUnique->nrow(); aux++) {
        mnAFull.push_back(SharedMatrix(new Matrix(basis_->nbf(), basis_->nbf())));
        mnAFull[aux]->set(AmnUnique->const_pointer()[aux]);
    }
    return mnAFull;
}

SharedMatrix MatPsi2::DFJK_InverseJHalfMetric() {
    DFJKException("DFJK_InverseJHalfMetric");
    return boost::static_pointer_cast<DFJK>(jk_)->GetInvJHalf();
}

void MatPsi2::RHF_Reset() {
    rhf_ = boost::shared_ptr<scf::RHF>(new scf::RHF(process_environment_, basis_));
    process_environment_.set_wavefunction(rhf_);
    rhf_->extern_finalize();
    rhf_.reset();
}

void MatPsi2::RHF_EnableMOM(int mom_start) {
    process_environment_.options.set_global_int("MOM_START", mom_start);
}

void MatPsi2::RHF_EnableDamping(double damping_percentage) {
    process_environment_.options.set_global_double("DAMPING_PERCENTAGE", damping_percentage);
}

void MatPsi2::RHF_DisableDIIS() {
    process_environment_.options.set_global_bool("DIIS", false);
    process_environment_.options.set_global_int("MAXITER", 500);
}

void MatPsi2::RHF_EnableDIIS() {
    process_environment_.options.set_global_int("DIIS", true);
    process_environment_.options.set_global_int("MAXITER", 100);
}

void MatPsi2::RHF_GuessSAD() {
    process_environment_.options.set_global_str("GUESS", "SAD");
}

void MatPsi2::RHF_GuessCore() {
    process_environment_.options.set_global_str("GUESS", "Core");
}

double MatPsi2::RHF_DoSCF() {
    if(Molecule_NumElectrons() % 2)
        throw PSIEXCEPTION("RHF_DoSCF: RHF can handle singlets only.");
    if(jk_ == NULL)
        JK_Initialize("PKJK");
    rhf_ = boost::shared_ptr<scf::RHF>(new scf::RHF(process_environment_, jk_));
    process_environment_.set_wavefunction(rhf_);
    return rhf_->compute_energy();
}

double MatPsi2::RHF_TotalEnergy() { 
    if(rhf_ == NULL) {
        throw PSIEXCEPTION("RHF_TotalEnergy: Hartree-Fock calculation has not been done.");
    }
    return rhf_->EHF(); 
}

SharedMatrix MatPsi2::RHF_Orbital() { 
    if(rhf_ == NULL) {
        throw PSIEXCEPTION("RHF_Orbital: Hartree-Fock calculation has not been done.");
    }
    return rhf_->Ca(); 
}

SharedVector MatPsi2::RHF_OrbitalEnergies() { 
    if(rhf_ == NULL) {
        throw PSIEXCEPTION("RHF_OrbitalEnergies: Hartree-Fock calculation has not been done.");
    }
    return rhf_->epsilon_a(); 
}

SharedMatrix MatPsi2::RHF_Density() { 
    if(rhf_ == NULL) {
        throw PSIEXCEPTION("RHF_Density: Hartree-Fock calculation has not been done.");
    }
    return rhf_->Da(); 
}

SharedMatrix MatPsi2::RHF_CoreHamiltonian() { // RHF_H leads to naming issues due to "ifdef RHF_H" or something 
    if(rhf_ == NULL) {
        throw PSIEXCEPTION("RHF_CoreHamiltonian: Hartree-Fock calculation has not been done.");
    }
    return rhf_->H(); 
}

SharedMatrix MatPsi2::RHF_JMatrix() { 
    if(rhf_ == NULL) {
        throw PSIEXCEPTION("RHF_JMatrix: Hartree-Fock calculation has not been done.");
    }
    return rhf_->J(); 
}

SharedMatrix MatPsi2::RHF_KMatrix() { 
    if(rhf_ == NULL) {
        throw PSIEXCEPTION("RHF_KMatrix: Hartree-Fock calculation has not been done.");
    }
    return rhf_->K(); 
}

SharedMatrix MatPsi2::RHF_FockMatrix() { 
    if(rhf_ == NULL) {
        throw PSIEXCEPTION("RHF_FockMatrix: Hartree-Fock calculation has not been done.");
    }
    return rhf_->Fa(); 
}

