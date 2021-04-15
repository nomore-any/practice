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

//!	�ڵ���
class CNode
{
public:

//!	ÿ���ڵ�����ɶȸ���������һ����̬��������������������κ�һ����������һ������
/*!	For 3D bar and solid elements, NDF = 3. For 3D beam or shell elements, NDF = 5 or 6 */
	const static unsigned int NDF = 3;

//!	����ڵ���
	unsigned int NodeNumber;

//!	����ڵ��XYZ��������
	double XYZ[3];

//!	���峤��Ϊ'�ڵ����ɶ���'������bcode�����ڴ洢Լ�������
//  ��CDomain���CalculateEquationNumber()ִ�й��󣬱�Ϊ�洢ȫ�ֽڵ���
/*!		0: The corresponding degree of freedom is active (defined in the global system) */
/*!		1: The corresponding degree of freedom in nonactive (not defined) */

	unsigned int bcode[NDF];

//!	���캯�������󱻹���ʱ���Զ�����X��Y��Z�����ֲ�����
	CNode(double X = 0, double Y = 0, double Z = 0);

//!	����һ��ifstream��������������ȡ��Ϣ
	bool Read(ifstream& Input);

//!	����COutputter��������������ӡ�ڵ���Ϣ
	void Write(COutputter& output);

//!	����COutputter�����������������ýڵ�3�����ɶȶ�Ӧ��ȫ�ֱ��
	void WriteEquationNo(COutputter& OutputFile);

//!	����һ��COutputter�����������������ȫ��λ�����飬����ýڵ��λ��
	void WriteNodalDisplacement(COutputter& OutputFile, double* Displacement);
};
