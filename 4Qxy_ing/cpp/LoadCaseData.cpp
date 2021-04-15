/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "LoadCaseData.h"

#include <iomanip>
#include <iostream>

using namespace std;

//���������Ȳ���
CLoadCaseData :: ~CLoadCaseData()
{
	delete [] node;
	delete [] dof;
	delete [] load;
}

//�����غ��������ڶ��������Ӧ����Ԫ�ظ�����node,dof,load����
void CLoadCaseData :: Allocate(unsigned int num)
{
	nloads = num;
	node = new unsigned int[nloads];
	dof = new unsigned int[nloads];
	load = new double[nloads];
}; 

//	��ȡnload�������ڴ�ռ䣬��ȡnode��dof��load������
bool CLoadCaseData :: Read(ifstream& Input)
{
	
	unsigned int NL;

	Input >> NL;

	Allocate(NL);

	for (unsigned int i = 0; i < NL; i++)
		Input >> node[i] >> dof[i] >> load[i];

	return true;
}

//	����غɲ���
void CLoadCaseData::Write(COutputter& output)
{
	for (unsigned int i = 0; i < nloads; i++)
		output << setw(7) << node[i] << setw(13) << dof[i]  << setw(19) << load[i] << endl;
}
