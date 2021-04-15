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

//! �غɹ��������࣬ÿ����һ���غɹ���������һ������Ķ���
class CLoadCaseData
{
public:

	unsigned int nloads;	//!�ù����ļ����غ���
	unsigned int* node;		//!< ��1~nloads,ÿ���غɶ��������ĸ��ڵ���
	unsigned int* dof;		//!< ��1~nloads,ÿ���غɶ��������ĸ����ɶ���
	double* load;			//!< ��1~nloads,ÿ���غɶ��ж��

public:
	//���캯�������ֶ�����ʼ��
	CLoadCaseData() : nloads(0), node(NULL), dof(NULL), load(NULL) {};
	~CLoadCaseData();

//!	���뼯���غ������ڶ�������ռ�
	void Allocate(unsigned int num);

//!	��ȡ�غ�����
	bool Read(ifstream& Input);

//!	����غ�����
	void Write(COutputter& output);
};
