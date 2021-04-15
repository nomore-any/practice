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
    XYZ[0] = X;		// 注意到X=Y=Z=0，相当于把坐标初始化为(0,0,0)
    XYZ[1] = Y;
    XYZ[2] = Z;
    
    bcode[0] = 0;	// 把约束情况全部初始化为0
    bcode[1] = 0;
    bcode[2] = 0;
};

//	以引用的方式输入一个ifstream对象，依靠这个对象提供的接口从文件中读取节点编号，节点约束情况，节点坐标。
//  之所以返回一个bool变量，是为了在CDomain的ReadNodalPoints()函数中便于判断输入格式是否正确。
bool CNode::Read(ifstream& Input)
{
	Input >> NodeNumber;	// node number
	Input >> bcode[0] >> bcode[1] >> bcode[2]
		  >> XYZ[0] >> XYZ[1] >> XYZ[2];

	return true;
}

//	以引用方式传入一个COutputter对象，以提供被写入文件的接口。
//  这个对象重载的运算符<<在COutputter中定义，一方面直接打印各种数据，另一方面把数据写进目标文件里。
void CNode::Write(COutputter& output)
{
	output << setw(9) << NodeNumber << setw(5) << bcode[0] << setw(5) << bcode[1] << setw(5) << bcode[2]
		   << setw(18) << XYZ[0] << setw(15) << XYZ[1] << setw(15) << XYZ[2] << endl;
}

//	输入一个COutputter对象以提供被写入文件的接口，输出该节点的三个自由度所对应的全局编号
void CNode::WriteEquationNo(COutputter& output)
{
	output << setw(9) << NodeNumber << "       ";

	for (unsigned int dof = 0; dof < CNode::NDF; dof++)	// Loop over for DOFs of node np
	{
		output << setw(5) << bcode[dof];
	}

	output << endl;
}

//	输入全局位移，输出节点位移
void CNode::WriteNodalDisplacement(COutputter& output, double* Displacement)
{
	output << setw(5) << NodeNumber << "        ";
	//对该节点的NDF个自由度进行循环
	//如果节点是被约束的，输出位移为0；如果节点不是被约束的，则根据bcode从全局位移中找到该节点对应的位移并输出
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
