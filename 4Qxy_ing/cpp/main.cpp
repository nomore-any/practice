/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Domain.h"
#include "Bar.h"
#include "Outputter.h"
#include "Clock.h"

using namespace std;
int main(int argc, char *argv[])
{
	//下面一大段都不用管
	if (argc != 2) //  Print help message
	{
	    cout << "Usage: stap++ InputFileName\n";
		exit(1);
	}

	string filename(argv[1]);
    size_t found = filename.find_last_of('.');

    // If the input file name is provided with an extension
    if (found != std::string::npos) {
        if (filename.substr(found) == ".dat")
            filename = filename.substr(0, found);
        else {
            // The input file name must has an extension of 'dat'
            cout << "*** Error *** Invalid file extension: "
                 << filename.substr(found+1) << endl;
            exit(1);
        }
    }

    string InFile = filename + ".dat";
	string OutFile = filename + ".out";

	 Clock timer;
    timer.Start();
	//上面一大段都不用管

	//获取指向CDomain的指针
	CDomain* FEMData = CDomain::GetInstance();

//  一条龙式地读取数据，不用管
	if (!FEMData->ReadData(InFile, OutFile))
	{
		cerr << "*** Error *** Data input failed!" << endl;
		exit(1);
	}
    
    double time_input = timer.ElapsedTime();

	//获取Outputter的指针，以便后面输出
    COutputter* Output = COutputter::GetInstance();

    if (!FEMData->GetMODEX())
    {
        *Output << "Data check completed !" << endl << endl;
        return 0;
    }

//  分配内存空间，计算与刚度阵相关的一些参数(比如列高)
	FEMData->AllocateMatrices();
    
//  组装并存储全局刚度阵
	FEMData->AssembleStiffnessMatrix();
    
    double time_assemble = timer.ElapsedTime();

//  把全局刚度阵传入一个CLDLTSolver类型的对象里以便后续求解
	CLDLTSolver* Solver = new CLDLTSolver(FEMData->GetStiffnessMatrix());
    
//  使用CLDLTSolver类里面的LDLT函数进行刚度阵的分解
    Solver->LDLT();

#ifdef _DEBUG_
    Output->PrintStiffnessMatrix();
#endif
        
//  对载荷工况的编号循环，每种工况都输出一个独立的位移和受力情况
    for (unsigned int lcase = 0; lcase < FEMData->GetNLCASE(); lcase++)
    {
//      Assemble righ-hand-side vector (force vector)
        FEMData->AssembleForce(lcase + 1);
            
//      Reduce right-hand-side force vector and back substitute
        Solver->BackSubstitution(FEMData->GetForce());

        *Output << " LOAD CASE" << setw(5) << lcase + 1 << endl << endl << endl;

#ifdef _DEBUG_
        Output->PrintDisplacement();
#endif
            //输出位移
        Output->OutputNodalDisplacement();

     //计算并输出应力和力
        Output->OutputElementStress();
    }

    double time_solution = timer.ElapsedTime();
    
    timer.Stop();
    
    *Output << "\n S O L U T I O N   T I M E   L O G   I N   S E C \n\n"
            << "     TIME FOR INPUT PHASE = " << time_input << endl
            << "     TIME FOR CALCULATION OF STIFFNESS MATRIX = " << time_assemble - time_input << endl
            << "     TIME FOR FACTORIZATION AND LOAD CASE SOLUTIONS = " << time_solution - time_assemble << endl << endl
            << "     T O T A L   S O L U T I O N   T I M E = " << time_solution << endl << endl;
	system("pause");
	return 0;
}
