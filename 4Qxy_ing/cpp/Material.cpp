/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Material.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

//	��ȡ���ϲ���
bool CBarMaterial::Read(ifstream& Input)
{
	Input >> nset;	//�������ϱ��(ÿ����Ԫ�鶼���Լ���1��ʼ�Ĳ��ϱ��)

	Input >> E >> Area;	// ��ȡE��A

	return true;
}

//	������ϲ���
void CBarMaterial::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << Area << endl;
}


//	��ȡ���ϲ���
bool CQ4Material::Read(ifstream& Input)
{
	Input >> nset;	//�������ϱ��(ÿ����Ԫ�鶼���Լ���1��ʼ�Ĳ��ϱ��)

	Input >> E >> Nu;	// ��ȡE��A�Ͳ��ɱ�

	return true;
}

//	������ϲ���
void CQ4Material::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << Nu << endl;
}
