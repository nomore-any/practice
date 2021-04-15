/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include "Outputter.h"

using namespace std;

//!	节点类
class CNode
{
public:

//!	每个节点的自由度个数，这是一个静态变量，属于类而不属于任何一个对象，且是一个常量
/*!	For 3D bar and solid elements, NDF = 3. For 3D beam or shell elements, NDF = 5 or 6 */
	const static unsigned int NDF = 3;

//!	定义节点编号
	unsigned int NodeNumber;

//!	定义节点的XYZ坐标数组
	double XYZ[3];

//!	定义长度为'节点自由度数'的数组bcode，用于存储约束情况；
//  在CDomain类的CalculateEquationNumber()执行过后，变为存储全局节点编号
/*!		0: The corresponding degree of freedom is active (defined in the global system) */
/*!		1: The corresponding degree of freedom in nonactive (not defined) */

	unsigned int bcode[NDF];

//!	构造函数，对象被构造时，自动传入X、Y、Z三个局部变量
	CNode(double X = 0, double Y = 0, double Z = 0);

//!	传入一个ifstream对象，用这个对象读取信息
	bool Read(ifstream& Input);

//!	传入COutputter对象，用这个对象打印节点信息
	void Write(COutputter& output);

//!	传入COutputter对象，用这个对象输出该节点3个自由度对应的全局编号
	void WriteEquationNo(COutputter& OutputFile);

//!	传入一个COutputter对象用于输出，传入全局位移数组，输出该节点的位移
	void WriteNodalDisplacement(COutputter& OutputFile, double* Displacement);
};
