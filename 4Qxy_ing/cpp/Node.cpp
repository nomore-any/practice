/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include <iostream>
#include <iomanip>

#include "Node.h"

CNode::CNode(double X, double Y, double Z)
{
    XYZ[0] = X;		// ע�⵽X=Y=Z=0���൱�ڰ������ʼ��Ϊ(0,0,0)
    XYZ[1] = Y;
    XYZ[2] = Z;
    
    bcode[0] = 0;	// ��Լ�����ȫ����ʼ��Ϊ0
    bcode[1] = 0;
    bcode[2] = 0;
};

//	�����õķ�ʽ����һ��ifstream����������������ṩ�Ľӿڴ��ļ��ж�ȡ�ڵ��ţ��ڵ�Լ��������ڵ����ꡣ
//  ֮���Է���һ��bool��������Ϊ����CDomain��ReadNodalPoints()�����б����ж������ʽ�Ƿ���ȷ��
bool CNode::Read(ifstream& Input)
{
	Input >> NodeNumber;	// node number
	Input >> bcode[0] >> bcode[1] >> bcode[2]
		  >> XYZ[0] >> XYZ[1] >> XYZ[2];

	return true;
}

//	�����÷�ʽ����һ��COutputter�������ṩ��д���ļ��Ľӿڡ�
//  ����������ص������<<��COutputter�ж��壬һ����ֱ�Ӵ�ӡ�������ݣ���һ���������д��Ŀ���ļ��
void CNode::Write(COutputter& output)
{
	output << setw(9) << NodeNumber << setw(5) << bcode[0] << setw(5) << bcode[1] << setw(5) << bcode[2]
		   << setw(18) << XYZ[0] << setw(15) << XYZ[1] << setw(15) << XYZ[2] << endl;
}

//	����һ��COutputter�������ṩ��д���ļ��Ľӿڣ�����ýڵ���������ɶ�����Ӧ��ȫ�ֱ��
void CNode::WriteEquationNo(COutputter& output)
{
	output << setw(9) << NodeNumber << "       ";

	for (unsigned int dof = 0; dof < CNode::NDF; dof++)	// Loop over for DOFs of node np
	{
		output << setw(5) << bcode[dof];
	}

	output << endl;
}

//	����ȫ��λ�ƣ�����ڵ�λ��
void CNode::WriteNodalDisplacement(COutputter& output, double* Displacement)
{
	output << setw(5) << NodeNumber << "        ";
	//�Ըýڵ��NDF�����ɶȽ���ѭ��
	//����ڵ��Ǳ�Լ���ģ����λ��Ϊ0������ڵ㲻�Ǳ�Լ���ģ������bcode��ȫ��λ�����ҵ��ýڵ��Ӧ��λ�Ʋ����
	for (unsigned int j = 0; j < NDF; j++)
	{
		if (bcode[j] == 0)
		{
			output << setw(18) << 0.0;
		}
		else
		{
			output << setw(18) << Displacement[bcode[j] - 1];
		}
	}

	output << endl;
}
