#include "mex.h"
#include "class_handle.hpp"
#include "MatPsi2.h"

using namespace std;
using namespace psi;
using namespace boost;

SharedMatrix InputMatrix(const mxArray*& Mat_m) {
    int nrow = mxGetM(Mat_m);
    int ncol = mxGetN(Mat_m);
    SharedMatrix Mat_c(new Matrix(nrow, ncol));
    double* Mat_m_pt = mxGetPr(Mat_m);
    double** Mat_c_pt = Mat_c->pointer();
    for(int i = 0; i < ncol; i++)
        for(int j = 0; j < nrow; j++)
            Mat_c_pt[j][i] = *Mat_m_pt++; // Matlab loops over a column first, but C++ loops over a row first 
    return Mat_c;
}

double InputScalar(const mxArray*& Mat_m) {
    double* Mat_m_pt = mxGetPr(Mat_m);
    return *Mat_m_pt;
}

void OutputMatrix(mxArray*& Mat_m, SharedMatrix Mat_c) {
    int nrow = Mat_c->nrow();
    int ncol = Mat_c->ncol();
    Mat_m = mxCreateDoubleMatrix( nrow, ncol, mxREAL);
    double* Mat_m_pt = mxGetPr(Mat_m);
    double** Mat_c_pt = Mat_c->pointer();
    for(int i = 0; i < ncol; i++)
        for(int j = 0; j < nrow; j++)
            *Mat_m_pt++ = Mat_c_pt[j][i];
}

void OutputVector(mxArray*& Mat_m, SharedVector Vec_c) {
    int dim = Vec_c->dim();
    Mat_m = mxCreateDoubleMatrix( 1, dim, mxREAL);
    double* Mat_m_pt = mxGetPr(Mat_m);
    double* Vec_c_pt = Vec_c->pointer();
    for(int i = 0; i < dim; i++) {
        *Mat_m_pt++ = *Vec_c_pt++;
    }
}

void OutputScalar(mxArray*& Mat_m, double scalar) {
    Mat_m = mxCreateDoubleMatrix( 1, 1, mxREAL);
    double* Mat_m_pt = mxGetPr(Mat_m);
    *Mat_m_pt = scalar;
}

void OutputVectorOfSymmMatrices(mxArray*& Mat_m, std::vector<SharedMatrix> vecOfMats) {
	int ndim1 = vecOfMats[0]->ncol();
	int ndim2 = vecOfMats[0]->nrow();
	int ndim3 = vecOfMats.size();
	mwSize dims[3] = {ndim1, ndim2, ndim3};
	Mat_m = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
	double* Mat_m_pt = mxGetPr(Mat_m);
	for(int idim3 = 0; idim3 < ndim3; idim3++) {
		double* tmp_pt = vecOfMats[idim3]->get_pointer();
		for(int i = 0; i < ndim1 * ndim2; i++) {
			Mat_m_pt[idim3*ndim1*ndim2 + i] = tmp_pt[i];
		}
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Get the command string
    char cmd[64];
    if (nrhs < 1 || mxGetString(prhs[0], cmd, sizeof(cmd)))
        mexErrMsgTxt("First input should be a command string less than 64 characters long.");
    
    // Constructor
    if (!strcmp("new", cmd)) {
        // Check parameters
        if (nlhs != 1)
            mexErrMsgTxt("MatPsi2 Constructor: One output expected.");
        if ( nrhs!=6 || !mxIsDouble(prhs[1]) || mxGetN(prhs[1]) != 4 || !mxIsChar(prhs[2]) || !mxIsDouble(prhs[3]) ||!mxIsDouble(prhs[4]) || !mxIsChar(prhs[5]) )
            mexErrMsgTxt("MatPsi2 Constructor: MatPsi(cartesian, basis_name, charge, multiplicity) input expected.");
        // Return a handle to a new C++ instance
        plhs[0] = convertPtr2Mat<MatPsi2>(new MatPsi2(
            InputMatrix(prhs[1]), 
            (std::string)mxArrayToString(prhs[2]), 
            (int)InputScalar(prhs[3]), 
            (int)InputScalar(prhs[4]),
            (std::string)mxArrayToString(prhs[5]) + "/"));
        return;
    }
    
    // Check there is a second input, which should be the class instance handle
    if (nrhs < 2)
        mexErrMsgTxt("Second input should be a class instance handle.");
    
    // Delete
    if (!strcmp("delete", cmd)) {
        // Destroy the C++ object
        destroyObject<MatPsi2>(prhs[1]);
        // Warn if other commands were ignored
        if (nlhs != 0 || nrhs != 2)
            mexWarnMsgTxt("Delete: Unexpected arguments ignored.");
        return;
    }
    
    // Get the class instance pointer from the second input
    MatPsi2* MatPsi_obj = convertMat2Ptr<MatPsi2>(prhs[1]);
    
    int nbf = MatPsi_obj->BasisSet_NumFunctions();
    
    //*** Call the various class methods 
    //*** Settings 
    if (!strcmp("Settings_MaxNumCPUCores", cmd)) {
        OutputScalar(plhs[0], (double)MatPsi_obj->Settings_MaxNumCPUCores());
        return;
    }
    if (!strcmp("Settings_MaxMemoryInGB", cmd)) {
        OutputScalar(plhs[0], MatPsi_obj->Settings_MaxMemoryInGB());
        return;
    }
    if (!strcmp("Settings_PsiDataDir", cmd)) {
        plhs[0] = mxCreateString((MatPsi_obj->Settings_PsiDataDir()).c_str());
        return;
    }
    if (!strcmp("Settings_TempDir", cmd)) {
        plhs[0] = mxCreateString((MatPsi_obj->Settings_TempDir()).c_str());
        return;
    }
    if (!strcmp("Settings_SetMaxNumCPUCores", cmd)) {
        if (nrhs == 2) {
            MatPsi_obj->Settings_SetMaxNumCPUCores(1);
            return;
        }
        if (nrhs!=3 ||  mxGetM(prhs[2])!=1 || mxGetN(prhs[2])!=1)
            mexErrMsgTxt("Settings_SetMaxNumCPUCores(ncores): Integer input expected.");
        MatPsi_obj->Settings_SetMaxNumCPUCores((int)InputScalar(prhs[2]));
        return;
    }
    if (!strcmp("Settings_SetMaxMemory", cmd)) {
        if (nrhs == 2) {
            MatPsi_obj->Settings_SetMaxMemory("1000mb");
            return;
        }
        if (nrhs!=3 || !mxIsChar(prhs[2]))
            mexErrMsgTxt("Settings_SetMaxMemory(\"memory\"): String input expected.");
        MatPsi_obj->Settings_SetMaxMemory((std::string)mxArrayToString(prhs[2]));
        return;
    }
    if (!strcmp("Settings_SetPsiDataDir", cmd)) {
        if (nrhs!=3 || !mxIsChar(prhs[2]))
            mexErrMsgTxt("Settings_SetPsiDataDir(\"psiDataDir\"): String input expected.");
        MatPsi_obj->Settings_SetPsiDataDir((std::string)mxArrayToString(prhs[2]));
        return;
    }
    
    //*** Molecule 
    if (!strcmp("Molecule_Fix", cmd)) {
        MatPsi_obj->Molecule_Fix();
        return;
    }
    if (!strcmp("Molecule_Free", cmd)) {
        MatPsi_obj->Molecule_Free();
        return;
    }
    if (!strcmp("Molecule_NumAtoms", cmd)) {
        OutputScalar(plhs[0], (double)MatPsi_obj->Molecule_NumAtoms());
        return;
    }
    if (!strcmp("Molecule_NumElectrons", cmd)) {
        OutputScalar(plhs[0], (double)MatPsi_obj->Molecule_NumElectrons());
        return;
    }
    if (!strcmp("Molecule_Geometry", cmd)) {
        OutputMatrix(plhs[0], MatPsi_obj->Molecule_Geometry());
        return;
    }
    if (!strcmp("Molecule_SetGeometry", cmd)) {
        // Check parameters
        if (nrhs!=3 || mxGetM(prhs[2]) != MatPsi_obj->Molecule_NumAtoms() || mxGetN(prhs[2]) != 3)
            mexErrMsgTxt("Molecule_SetGeometry(newGeom): NumAtoms by 3 matrix input expected.");
        // Call the method
        MatPsi_obj->Molecule_SetGeometry(InputMatrix(prhs[2]));
        return;
    }
    if (!strcmp("Molecule_AtomicNumbers", cmd)) {
        OutputVector(plhs[0], MatPsi_obj->Molecule_AtomicNumbers());
        return;
    }
    if (!strcmp("Molecule_NucRepEnergy", cmd)) {
        OutputScalar(plhs[0], MatPsi_obj->Molecule_NucRepEnergy());
        return;
    }
    if (!strcmp("Molecule_ChargeMult", cmd)) {
        OutputVector(plhs[0], MatPsi_obj->Molecule_ChargeMult());
        return;
    }
    
    //*** BasisSet 
    if (!strcmp("BasisSet_Name", cmd)) {
        plhs[0] = mxCreateString((MatPsi_obj->BasisSet_Name()).c_str());
        return;
    }
    if (!strcmp("BasisSet_IsSpherical", cmd)) {
        OutputScalar(plhs[0], (double)MatPsi_obj->BasisSet_IsSpherical());
        return;
    }
    if (!strcmp("BasisSet_NumFunctions", cmd)) {
        OutputScalar(plhs[0], (double)MatPsi_obj->BasisSet_NumFunctions());
        return;
    }
    if (!strcmp("BasisSet_NumShells", cmd)) {
        OutputScalar(plhs[0], (double)MatPsi_obj->BasisSet_NumShells());
        return;
    }
    if (!strcmp("BasisSet_ShellTypes", cmd)) {
        OutputVector(plhs[0], MatPsi_obj->BasisSet_ShellTypes());
        return;
    }
    if (!strcmp("BasisSet_ShellNumPrimitives", cmd)) {
        OutputVector(plhs[0], MatPsi_obj->BasisSet_ShellNumPrimitives());
        return;
    }
    if (!strcmp("BasisSet_ShellNumFunctions", cmd)) {
        OutputVector(plhs[0], MatPsi_obj->BasisSet_ShellNumFunctions());
        return;
    }
    if (!strcmp("BasisSet_ShellToCenter", cmd)) {
        SharedVector shell2centerVec = MatPsi_obj->BasisSet_ShellToCenter();
        for(int i = 0; i < shell2centerVec->dim(); i++)
            shell2centerVec->add(i, 1.0); // + 1 convert C++ convention to Matlab convention 
        OutputVector(plhs[0], shell2centerVec);
        return;
    }
    if (!strcmp("BasisSet_FuncToCenter", cmd)) {
        SharedVector func2centerVec = MatPsi_obj->BasisSet_FuncToCenter();
        for(int i = 0; i < func2centerVec->dim(); i++)
            func2centerVec->add(i, 1.0); // + 1 convert C++ convention to Matlab convention 
        OutputVector(plhs[0], func2centerVec);
        return;
    }
    if (!strcmp("BasisSet_FuncToShell", cmd)) {
        SharedVector func2shellVec = MatPsi_obj->BasisSet_FuncToShell();
        for(int i = 0; i < func2shellVec->dim(); i++)
            func2shellVec->add(i, 1.0); // + 1 convert C++ convention to Matlab convention 
        OutputVector(plhs[0], func2shellVec);
        return;
    }
    if (!strcmp("BasisSet_FuncToAngular", cmd)) {
        OutputVector(plhs[0], MatPsi_obj->BasisSet_FuncToAngular());
        return;
    }
    if (!strcmp("BasisSet_PrimExp", cmd)) {
        OutputVector(plhs[0], MatPsi_obj->BasisSet_PrimExp());
        return;
    }
    if (!strcmp("BasisSet_PrimCoeffUnnorm", cmd)) {
        OutputVector(plhs[0], MatPsi_obj->BasisSet_PrimCoeffUnnorm());
        return;
    }
    
    
    //*** Integral 
    if (!strcmp("Integrals_Overlap", cmd)) {
        OutputMatrix(plhs[0], MatPsi_obj->Integrals_Overlap());
        return;
    }
    if (!strcmp("Integrals_Kinetic", cmd)) {
        OutputMatrix(plhs[0], MatPsi_obj->Integrals_Kinetic());
        return;
    }
    if (!strcmp("Integrals_Potential", cmd)) {
        OutputMatrix(plhs[0], MatPsi_obj->Integrals_Potential());
        return;
    }
    if (!strcmp("Integrals_PotentialEachCore", cmd)) {
        OutputVectorOfSymmMatrices(plhs[0], MatPsi_obj->Integrals_PotentialEachCore());
        return;
    }
    if (!strcmp("Integrals_PotentialPtQ", cmd)) {
        // Check parameters
        if (nrhs!=3)
            mexErrMsgTxt("Integrals_PotentialPtQ(Zxyz_mat): (number of point charges) by 4 matrix input expected.");
        if (mxGetN(prhs[2]) != 4)
            mexErrMsgTxt("Integrals_PotentialPtQ: Zxyz list matrix dimension does not agree.");
        // Call the method
        OutputMatrix(plhs[0], MatPsi_obj->Integrals_PotentialPtQ(InputMatrix(prhs[2])));
        return;
    }
    if (!strcmp("Integrals_Dipole", cmd)) {
        std::vector<SharedMatrix> dipole = MatPsi_obj->Integrals_Dipole();
        if(nlhs == 3) {
            OutputMatrix(plhs[0], dipole[0]);
            OutputMatrix(plhs[1], dipole[1]);
            OutputMatrix(plhs[2], dipole[2]);
            return;
        }
        OutputVectorOfSymmMatrices(plhs[0], dipole);
        return;
    }
    if (!strcmp("Integrals_ijkl", cmd)) {
        // Check parameters
        if (nrhs!=6 || !mxIsDouble(prhs[2]) || !mxIsDouble(prhs[3]) || !mxIsDouble(prhs[4]) || !mxIsDouble(prhs[5]))
            mexErrMsgTxt("Integrals_ijkl(i, j, k, l): 4 integers input expected.");
        // Call the method
        int ind[4];
        for(int i = 0; i < 4; i++) {
            ind[i] = (int)mxGetScalar(prhs[2+i]) - 1; // -1 convert Matlab convention to C++ convention
            if(ind[i] < 0 || ind[i] >= nbf)
                mexErrMsgTxt("Integrals_ijkl: Required index not within scale.");
        }
        OutputScalar(plhs[0], MatPsi_obj->Integrals_ijkl(ind[0], ind[1], ind[2], ind[3]));
        return;
    }
    if (!strcmp("Integrals_NumUniqueTEIs", cmd)) {
        OutputScalar(plhs[0], (double)MatPsi_obj->Integrals_NumUniqueTEIs());
        return;
    }
    if (!strcmp("Integrals_AllUniqueTEIs", cmd)) {
        plhs[0] = mxCreateDoubleMatrix( 1, MatPsi_obj->Integrals_NumUniqueTEIs(), mxREAL);
        double* matpt = mxGetPr(plhs[0]);
        MatPsi_obj->Integrals_AllUniqueTEIs(matpt);
        return;
    }
    if (!strcmp("Integrals_AllTEIs", cmd)) {
        mwSize dims[4] = {nbf, nbf, nbf, nbf};
        plhs[0] = mxCreateNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL);
        double* matpt = mxGetPr(plhs[0]);
        MatPsi_obj->Integrals_AllTEIs(matpt);
        return;
    }
    if (!strcmp("Integrals_IndicesForK", cmd)) {
        plhs[0] = mxCreateDoubleMatrix( 1, MatPsi_obj->Integrals_NumUniqueTEIs(), mxREAL);
        double* matpt1 = mxGetPr(plhs[0]);
        plhs[1] = mxCreateDoubleMatrix( 1, MatPsi_obj->Integrals_NumUniqueTEIs(), mxREAL);
        double* matpt2 = mxGetPr(plhs[1]);
        MatPsi_obj->Integrals_IndicesForK(matpt1, matpt2);
        return;
    }
    
    //*** JK related  
    if (!strcmp("JK_Initialize", cmd)) {
        if ( (nrhs!=3 && nrhs!=4) || !mxIsChar(prhs[2]) )
            mexErrMsgTxt("JK_Initialize(\"jktype\"): String input expected.");
        if (nrhs==3) {
            MatPsi_obj->JK_Initialize((std::string)mxArrayToString(prhs[2]));
        }else {
            MatPsi_obj->JK_Initialize((std::string)mxArrayToString(prhs[2]), (std::string)mxArrayToString(prhs[3]));
        }
        return;
    }
    if (!strcmp("JK_Type", cmd)) {
        plhs[0] = mxCreateString((MatPsi_obj->JK_Type()).c_str());
        return;
    }    
    if (!strcmp("JK_DensToJ", cmd)) {
        std::vector<SharedMatrix> vecOfJMats;
        if (nrhs==3 && mxGetM(prhs[2]) == nbf && mxGetN(prhs[2]) == nbf)
            vecOfJMats = MatPsi_obj->JK_DensToJ(InputMatrix(prhs[2]));
        else if (nrhs==4 && mxGetM(prhs[2]) == nbf && mxGetN(prhs[2]) == nbf && mxGetM(prhs[3]) == nbf && mxGetN(prhs[3]) == nbf)
            vecOfJMats = MatPsi_obj->JK_DensToJ(InputMatrix(prhs[2]), InputMatrix(prhs[3]));
        else
            mexErrMsgTxt("JK_DensToJ(densAlpha, densBeta): 1 or 2 nbf by nbf matrix(ces) input expected.");
        OutputVectorOfSymmMatrices(plhs[0], vecOfJMats);
        return;
    }
    if (!strcmp("JK_DensToK", cmd)) {
        std::vector<SharedMatrix> vecOfKMats;
        if (nrhs==3 && mxGetM(prhs[2]) == nbf && mxGetN(prhs[2]) == nbf)
            vecOfKMats = MatPsi_obj->JK_DensToK(InputMatrix(prhs[2]));
        else if (nrhs==4 && mxGetM(prhs[2]) == nbf && mxGetN(prhs[2]) == nbf && mxGetM(prhs[3]) == nbf && mxGetN(prhs[3]) == nbf)
            vecOfKMats = MatPsi_obj->JK_DensToK(InputMatrix(prhs[2]), InputMatrix(prhs[3]));
        else
            mexErrMsgTxt("JK_DensToK(densAlpha, densBeta): 1 or 2 nbf by nbf matrix(ces) input expected.");
        OutputVectorOfSymmMatrices(plhs[0], vecOfKMats);
        return;
    }
    if (!strcmp("JK_OccOrbToJ", cmd)) {
        std::vector<SharedMatrix> vecOfJMats;
        if (nrhs==3 && mxGetM(prhs[2]) == nbf)
            vecOfJMats = MatPsi_obj->JK_OccOrbToJ(InputMatrix(prhs[2]));
        else if (nrhs==4 && mxGetM(prhs[2]) == nbf && mxGetM(prhs[3]) == nbf)
            vecOfJMats = MatPsi_obj->JK_OccOrbToJ(InputMatrix(prhs[2]), InputMatrix(prhs[3]));
        else
            mexErrMsgTxt("JK_OccOrbToJ(occOrbAlpha, occOrbBeta): 1 or 2 nbf by any matrix(ces) input expected.");
        OutputVectorOfSymmMatrices(plhs[0], vecOfJMats);
        return;
    }
    if (!strcmp("JK_OccOrbToK", cmd)) {
        std::vector<SharedMatrix> vecOfKMats;
        if (nrhs==3 && mxGetM(prhs[2]) == nbf)
            vecOfKMats = MatPsi_obj->JK_OccOrbToK(InputMatrix(prhs[2]));
        else if (nrhs==4 && mxGetM(prhs[2]) == nbf && mxGetM(prhs[3]) == nbf)
            vecOfKMats = MatPsi_obj->JK_OccOrbToK(InputMatrix(prhs[2]), InputMatrix(prhs[3]));
        else
            mexErrMsgTxt("JK_OccOrbToK(occOrbAlpha, occOrbBeta): 1 or 2 nbf by any matrix(ces) input expected.");
        OutputVectorOfSymmMatrices(plhs[0], vecOfKMats);
        return;
    }
    if (!strcmp("JK_CalcAllFromDens", cmd)) {
        if (nrhs==3 && mxGetM(prhs[2]) == nbf && mxGetN(prhs[2]) == nbf)
            MatPsi_obj->JK_CalcAllFromDens(InputMatrix(prhs[2]));
        else if (nrhs==4 && mxGetM(prhs[2]) == nbf && mxGetN(prhs[2]) == nbf && mxGetM(prhs[3]) == nbf && mxGetN(prhs[3]) == nbf)
            MatPsi_obj->JK_CalcAllFromDens(InputMatrix(prhs[2]), InputMatrix(prhs[3]));
        else
            mexErrMsgTxt("JK_CalcAllFromDens(densAlpha, densBeta): 1 or 2 nbf by nbf matrix(ces) input expected.");
        return;
    }
    if (!strcmp("JK_CalcAllFromOccOrb", cmd)) {
        if (nrhs==3 && mxGetM(prhs[2]) == nbf)
            MatPsi_obj->JK_CalcAllFromOccOrb(InputMatrix(prhs[2]));
        else if (nrhs==4 && mxGetM(prhs[2]) == nbf && mxGetM(prhs[3]) == nbf)
            MatPsi_obj->JK_CalcAllFromOccOrb(InputMatrix(prhs[2]), InputMatrix(prhs[3]));
        else
            mexErrMsgTxt("JK_CalcAllFromOccOrb(occOrbAlpha, occOrbBeta): 1 or 2 nbf by any matrix(ces) input expected.");
        return;
    }
    if (!strcmp("JK_RetrieveJ", cmd)) {
        OutputVectorOfSymmMatrices(plhs[0], MatPsi_obj->JK_RetrieveJ());
        return;
    }
    if (!strcmp("JK_RetrieveK", cmd)) {
        OutputVectorOfSymmMatrices(plhs[0], MatPsi_obj->JK_RetrieveK());
        return;
    }
    if (!strcmp("JK_DFTensor_AuxPriPairs", cmd)) {
        // Call the method
        OutputMatrix(plhs[0], MatPsi_obj->JK_DFTensor_AuxPriPairs());
        return;
    }
    if (!strcmp("JK_DFTensor_AuxPriPri", cmd)) {
        OutputVectorOfSymmMatrices(plhs[0], MatPsi_obj->JK_DFTensor_AuxPriPri());
        return;
    }
    if (!strcmp("JK_DFMetric_InvJHalf", cmd)) {
        // Call the method
        OutputMatrix(plhs[0], MatPsi_obj->JK_DFMetric_InvJHalf());
        return;
    }
    
    //*** DFT related 
    if (!strcmp("DFT_Initialize", cmd)) {
        if ( nrhs!=3 || !mxIsChar(prhs[2]) )
            mexErrMsgTxt("DFT_Initialize(\"functionalName\"): String input expected.");
        MatPsi_obj->DFT_Initialize((std::string)mxArrayToString(prhs[2]));
        return;
    }
    if (!strcmp("DFT_DensToV", cmd)) {
        std::vector<SharedMatrix> dftPotArray;
        if (nrhs==3 && mxGetM(prhs[2]) == nbf)
            dftPotArray = MatPsi_obj->DFT_DensToV(InputMatrix(prhs[2]));
        else if (nrhs==4 && mxGetM(prhs[2]) == nbf && mxGetN(prhs[2]) == nbf && mxGetM(prhs[3]) == nbf && mxGetN(prhs[3]) == nbf)
            dftPotArray = MatPsi_obj->DFT_DensToV(InputMatrix(prhs[2]), InputMatrix(prhs[3]));
        else
            mexErrMsgTxt("DFT_DensToV(densAlpha, densBeta): 1 or 2 nbf by nbf matrix(ces) input expected.");
        OutputVectorOfSymmMatrices(plhs[0], dftPotArray);
        return;
    }
    if (!strcmp("DFT_OccOrbToV", cmd)) {
        std::vector<SharedMatrix> dftPotArray;
        if (nrhs==3 && mxGetM(prhs[2]) == nbf)
            dftPotArray = MatPsi_obj->DFT_OccOrbToV(InputMatrix(prhs[2]));
        else if (nrhs==4 && mxGetM(prhs[2]) == nbf && mxGetM(prhs[3]) == nbf)
            dftPotArray = MatPsi_obj->DFT_OccOrbToV(InputMatrix(prhs[2]), InputMatrix(prhs[3]));
        else
            mexErrMsgTxt("DFT_OccOrbToV(occOrbAlpha, occOrbBeta): 1 or 2 nbf by any matrix(ces) input expected.");
        OutputVectorOfSymmMatrices(plhs[0], dftPotArray);
        return;
    }
    if (!strcmp("DFT_EnergyXC", cmd)) {
        OutputScalar(plhs[0], MatPsi_obj->DFT_EnergyXC());
        return;
    }
    
    //*** SCF related 
    if (!strcmp("SCF_SetSCFType", cmd)) {
        if ( nrhs!=3 || !mxIsChar(prhs[2]) )
            mexErrMsgTxt("SCF_SetSCFType(\"scfType\"): String input expected.");
        MatPsi_obj->SCF_SetSCFType((std::string)mxArrayToString(prhs[2]));
        return;
    }
    if (!strcmp("SCF_SetGuessOrb", cmd)) {
        // Check parameters
        if (nrhs==3 && mxGetM(prhs[2]) == nbf)
            MatPsi_obj->SCF_SetGuessOrb(InputMatrix(prhs[2]));
        else if (nrhs==4 && mxGetM(prhs[2]) == nbf && mxGetM(prhs[3]) == nbf)
            MatPsi_obj->SCF_SetGuessOrb(InputMatrix(prhs[2]), InputMatrix(prhs[3]));
        else
            mexErrMsgTxt("SCF_SetGuessOrb(occOrbAlpha, occOrbBeta): 1 or 2 nbf by any matrix(ces) input expected.");
        return;
    }
    
    
    if (!strcmp("SCF_RunSCF", cmd)) {
        OutputScalar(plhs[0], MatPsi_obj->SCF_RunSCF());
        return;
    }
    if (!strcmp("SCF_EnableMOM", cmd)) {
        if (nrhs == 2) {
            MatPsi_obj->SCF_EnableMOM(20);
            return;
        }
        if (nrhs!=3 || mxGetM(prhs[2])!=1 || mxGetN(prhs[2])!=1)
            mexErrMsgTxt("SCF_EnableMOM(mom_start): Integer input expected.");
        MatPsi_obj->SCF_EnableMOM((int)InputScalar(prhs[2]));
        return;
    }
    if (!strcmp("SCF_DisableMOM", cmd)) {
        MatPsi_obj->SCF_EnableMOM(0);
        return;
    }
    if (!strcmp("SCF_EnableDamping", cmd)) {
        if (nrhs == 2) {
            MatPsi_obj->SCF_EnableDamping(20.0);
            return;
        }
        if (nrhs!=3 || mxGetM(prhs[2])!=1 || mxGetN(prhs[2])!=1)
            mexErrMsgTxt("SCF_EnableDamping(damping_percentage): 1 double input expected.");
        MatPsi_obj->SCF_EnableDamping(InputScalar(prhs[2]));
        return;
    }
    if (!strcmp("SCF_DisableDamping", cmd)) {
        MatPsi_obj->SCF_EnableDamping(0.0);
        return;
    }
    if (!strcmp("SCF_EnableDIIS", cmd)) {
        MatPsi_obj->SCF_EnableDIIS();
        return;
    }
    if (!strcmp("SCF_DisableDIIS", cmd)) {
        MatPsi_obj->SCF_DisableDIIS();
        return;
    }
    if (!strcmp("SCF_GuessSAD", cmd)) {
        MatPsi_obj->SCF_GuessSAD();
        return;
    }
    if (!strcmp("SCF_GuessCore", cmd)) {
        MatPsi_obj->SCF_GuessCore();
        return;
    }
    if (!strcmp("SCF_TotalEnergy", cmd)) {
        OutputScalar(plhs[0], MatPsi_obj->SCF_TotalEnergy());
        return;
    }
    if (!strcmp("SCF_OrbitalAlpha", cmd)) {
        OutputMatrix(plhs[0], MatPsi_obj->SCF_OrbitalAlpha());
        return;
    }
    if (!strcmp("SCF_OrbitalBeta", cmd)) {
        OutputMatrix(plhs[0], MatPsi_obj->SCF_OrbitalBeta());
        return;
    }
    if (!strcmp("SCF_OrbEigValAlpha", cmd)) {
        OutputVector(plhs[0], MatPsi_obj->SCF_OrbEigValAlpha());
        return;
    }
    if (!strcmp("SCF_OrbEigValBeta", cmd)) {
        OutputVector(plhs[0], MatPsi_obj->SCF_OrbEigValBeta());
        return;
    }
    if (!strcmp("SCF_DensityAlpha", cmd)) {
        OutputMatrix(plhs[0], MatPsi_obj->SCF_DensityAlpha());
        return;
    }
    if (!strcmp("SCF_DensityBeta", cmd)) {
        OutputMatrix(plhs[0], MatPsi_obj->SCF_DensityBeta());
        return;
    }
    if (!strcmp("SCF_FockAlpha", cmd)) {
        OutputMatrix(plhs[0], MatPsi_obj->SCF_FockAlpha());
        return;
    }
    if (!strcmp("SCF_FockBeta", cmd)) {
        OutputMatrix(plhs[0], MatPsi_obj->SCF_FockBeta());
        return;
    }
    if (!strcmp("SCF_Gradient", cmd)) {
        OutputMatrix(plhs[0], MatPsi_obj->SCF_Gradient());
        return;
    }
    if (!strcmp("SCF_GuessDensity", cmd)) {
        OutputMatrix(plhs[0], MatPsi_obj->SCF_GuessDensity());
        return;
    }
    if (!strcmp("SCF_RHF_J", cmd)) {
        OutputMatrix(plhs[0], MatPsi_obj->SCF_RHF_J());
        return;
    }
    if (!strcmp("SCF_RHF_K", cmd)) {
        OutputMatrix(plhs[0], MatPsi_obj->SCF_RHF_K());
        return;
    }
    
    // Got here, so command not recognized
    mexErrMsgTxt("Command not recognized.");
}

