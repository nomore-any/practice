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
	//����һ��ζ����ù�
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
	//����һ��ζ����ù�

	//��ȡָ��CDomain��ָ��
	CDomain* FEMData = CDomain::GetInstance();

//  һ����ʽ�ض�ȡ���ݣ����ù�
	if (!FEMData->ReadData(InFile, OutFile))
	{
		cerr << "*** Error *** Data input failed!" << endl;
		exit(1);
	}
    
    double time_input = timer.ElapsedTime();

	//��ȡOutputter��ָ�룬�Ա�������
    COutputter* Output = COutputter::GetInstance();

    if (!FEMData->GetMODEX())
    {
        *Output << "Data check completed !" << endl << endl;
        return 0;
    }

//  �����ڴ�ռ䣬������ն�����ص�һЩ����(�����и�)
	FEMData->AllocateMatrices();
    
//  ��װ���洢ȫ�ָն���
	FEMData->AssembleStiffnessMatrix();
    
    double time_assemble = timer.ElapsedTime();

//  ��ȫ�ָն�����һ��CLDLTSolver���͵Ķ������Ա�������
	CLDLTSolver* Solver = new CLDLTSolver(FEMData->GetStiffnessMatrix());
    
//  ʹ��CLDLTSolver�������LDLT�������иն���ķֽ�
    Solver->LDLT();

#ifdef _DEBUG_
    Output->PrintStiffnessMatrix();
#endif
        
//  ���غɹ����ı��ѭ����ÿ�ֹ��������һ��������λ�ƺ��������
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
            //���λ��
        Output->OutputNodalDisplacement();

     //���㲢���Ӧ������
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
