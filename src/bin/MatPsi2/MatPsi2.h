#include <libmints/mints.h>
#include <libfock/jk.h>
#include <libfock/v.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <boost/shared_array.hpp>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>
#include <libscf_solver/rhf.h>
#include <libscf_solver/ks.h>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include <libscfgrad/scf_grad.h>


using namespace std;
using namespace psi;
using namespace boost;

class MatPsi2 {
protected:

    Process::Environment process_environment_;
    boost::shared_ptr<LocalCommWrapper> worldcomm_;
    
    boost::shared_ptr<PSIO> psio_;
    
    boost::shared_ptr<Molecule> molecule_;
    boost::shared_ptr<BasisSet> basis_;
    boost::shared_ptr<IntegralFactory> intfac_;
    boost::shared_ptr<TwoBodyAOInt> eri_;
    boost::shared_ptr<MatrixFactory> matfac_;
    boost::shared_ptr<JK> jk_;
    boost::shared_ptr<VBase> dftPotential_;
    boost::shared_ptr<scf::HF> wfn_;
    
    SharedMatrix guessOrbital_;
    
    // create basis object 
    void create_basis();
    
    // create basis object and one & two electron integral factories 
    void create_basis_and_integral_factories();
    
    void create_wfn();
    
    // exception function for DFJK utilities
    void jk_DFException(std::string functionName);
    
public:
    // constructor
    MatPsi2(SharedMatrix cartesian, const std::string& basisname, 
        int charge, int multiplicity, const std::string& path);
    
    // destructor 
    virtual ~MatPsi2();
    
    // CPU and memory controll 
    int Settings_MaxNumCPUCores() { return process_environment_.get_n_threads(); }
    double Settings_MaxMemoryInGB() { return (double)process_environment_.get_memory() / 1E+9; }
    std::string Settings_PsiDataDir() { return process_environment_("PSIDATADIR"); }
    std::string Settings_TempDir() { return psio_->_psio_manager_->get_file_path(0); }
    void Settings_SetMaxNumCPUCores(int ncores);
    void Settings_SetMaxMemory(std::string);
    void Settings_SetPsiDataDir(std::string path) { process_environment_.set("PSIDATADIR", path); }
    
    
    
    //*** Molecule properties 
    int Molecule_NumAtoms() { return molecule_->natom(); } // number of atoms 
    int Molecule_NumElectrons(); // number of electrons 
    SharedMatrix Molecule_Geometry() { return molecule_->geometry().clone(); } // geometry in Bohr 
    void Molecule_SetGeometry(SharedMatrix newGeom); // set a new geometry in Bohr 
    double Molecule_NucRepEnergy() { return molecule_->nuclear_repulsion_energy(); } // nuclear repulsion energy 
    SharedVector Molecule_AtomicNumbers(); // atomic number list vector 
    void Molecule_SetChargeMult(int charge, int mult);
    SharedVector Molecule_ChargeMult();
    
    //*** Molecule operations 
    void Molecule_Fix();
    void Molecule_Free();
    
    
    //*** Basis set properties 
    const std::string& BasisSet_Name() { return basis_->name(); } // basis set name string 
    void BasisSet_SetBasisSet(const std::string& basisname); // set a new basis set 
    bool BasisSet_IsSpherical() { return basis_->has_puream(); }
    int BasisSet_NumFunctions() { return basis_->nbf(); } // number of basis functions 
    int BasisSet_NumShells() { return basis_->nshell(); }
    SharedVector BasisSet_ShellTypes();
    SharedVector BasisSet_ShellNumPrimitives();
    SharedVector BasisSet_ShellNumFunctions();
    SharedVector BasisSet_ShellToCenter();
    SharedVector BasisSet_FuncToCenter(); // map basis function to the index of the atom it is centred on 
    SharedVector BasisSet_FuncToShell();
    SharedVector BasisSet_FuncToAngular(); // map basis function number to its angular momentum 
    SharedVector BasisSet_PrimExp();
    SharedVector BasisSet_PrimCoeffUnnorm();
    
    
    //*** Integral package
    SharedMatrix Integrals_Overlap(); // overlap matrix S <i|j>
    SharedMatrix Integrals_Kinetic(); // kinetic energy matrix KE 
    SharedMatrix Integrals_Potential(); // total potential energy matrix EN <i|sum(1/R)|j>
    std::vector<SharedMatrix> Integrals_Dipole(); // dipole matrices <i|x|j>, <i|y|j>, <i|z|j>
    std::vector<SharedMatrix> Integrals_PotentialEachCore(); // atom-separated EN 
    SharedMatrix Integrals_PotentialPtQ(SharedMatrix Zxyz_list); // compute from a given point charge list the environment potential energy matrix ENVI 
    int Integrals_NumUniqueTEIs(); // number of unique TEIs 
    double Integrals_ijkl(int i, int j, int k, int l); // (ij|kl), chemist's notation 
    // ## HIGH MEMORY COST METHODS ## 
    void Integrals_AllUniqueTEIs(double*); // all unique TEIs in a vector 
    void Integrals_AllTEIs(double*); // all (repetitive) TEIs in a 4D-array 
    void Integrals_IndicesForK(double*, double*); // pre-arrange TEI vectors for forming K 
    // ## HIGH MEMORY COST METHODS ## 
    
    
    //*** JK related
    // use different types of JK 
    void JK_Initialize(std::string jktype, std::string auxiliaryBasisSetName = "CC-PVDZ-JKFIT");
    const std::string& JK_Type();
    
    // methods computing J/K 
    std::vector<SharedMatrix> JK_DensToJ(SharedMatrix, SharedMatrix = SharedMatrix());
    std::vector<SharedMatrix> JK_DensToK(SharedMatrix, SharedMatrix = SharedMatrix());
    std::vector<SharedMatrix> JK_OccOrbToJ(SharedMatrix, SharedMatrix = SharedMatrix());
    std::vector<SharedMatrix> JK_OccOrbToK(SharedMatrix, SharedMatrix = SharedMatrix());
    void JK_CalcAllFromDens(SharedMatrix, SharedMatrix = SharedMatrix());
    void JK_CalcAllFromOccOrb(SharedMatrix, SharedMatrix = SharedMatrix());
    std::vector<SharedMatrix> JK_RetrieveJ();
    std::vector<SharedMatrix> JK_RetrieveK();
    
    // specially for density-fitting JK
    SharedMatrix JK_DFTensor_AuxPriPairs();
    std::vector<SharedMatrix> JK_DFTensor_AuxPriPri();
    SharedMatrix JK_DFMetric_InvJHalf();
    
    
    void DFT_Initialize(std::string);
    std::vector<SharedMatrix> DFT_DensToV(SharedMatrix, SharedMatrix = SharedMatrix());
    std::vector<SharedMatrix> DFT_OccOrbToV(SharedMatrix, SharedMatrix = SharedMatrix());
    double DFT_EnergyXC();
    
    
    //*** SCF related
    // method of doing RHF calculations 
    void SCF_SetSCFType(std::string scfType);
    void SCF_SetGuessOrb(SharedMatrix guessOrb);
    double SCF_RunSCF();
    
    // methods controlling RHF algorithm 
    void SCF_EnableMOM(int mom_start);
    void SCF_EnableDamping(double damping_percentage);
    void SCF_EnableDIIS();
    void SCF_DisableDIIS();
    void SCF_GuessSAD();
    void SCF_GuessCore();
    
    // methods extracting restricted Hartree-Fock results
    double SCF_TotalEnergy();
    SharedMatrix SCF_OrbitalAlpha();
    SharedMatrix SCF_OrbitalBeta();
    SharedVector SCF_OrbEigValAlpha();
    SharedVector SCF_OrbEigValBeta();
    SharedMatrix SCF_DensityAlpha();
    SharedMatrix SCF_DensityBeta();
    SharedMatrix SCF_CoreHamiltonian();
    SharedMatrix SCF_FockAlpha();
    SharedMatrix SCF_FockBeta();
    SharedMatrix SCF_GuessDensity();
    
    SharedMatrix SCF_Gradient();
    
    SharedMatrix SCF_RHF_J();
    SharedMatrix SCF_RHF_K();
    
};
