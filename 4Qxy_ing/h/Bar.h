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

#include "Element.h"

using namespace std;

//! ��Element��̳е�Bar��
class CBar : public CElement
{
public:

//!	���캯��
	CBar();

//!	��������
	~CBar();

//!	��ȡBar��Ԫ������
	virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

//!	��ӡBar��Ԫ������
	virtual void Write(COutputter& output);

//!	���㵥Ԫ�ն��󣬴浽�����������
	virtual void ElementStiffness(double* Matrix);

//!	���㵥ԪӦ��
	virtual void ElementStress(double* stress, double* Displacement);
};
