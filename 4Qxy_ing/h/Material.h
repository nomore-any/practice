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
//����һ�����ϵĸ��࣬��������ж��������в��϶����е�����
//!	Material base class which only define one data member
/*!	All type of material classes should be derived from this base class */
class CMaterial
{
public:

	unsigned int nset;	//��������ϵı�ţ����ϱ��1,2,3...����ElementGroup�仯��ˢ�µ�
	
	double E;  //!������ϵ�����ģ��

public:

//! ����������
    virtual ~CMaterial() {};

//!	�麯�������ڶ�ȡ���ϲ���
	virtual bool Read(ifstream& Input) = 0;

//!	�麯��������������ϲ���
    virtual void Write(COutputter& output) = 0;

};

//!	Bar���͵�Ԫ��Ӧ�Ĳ�����
class CBarMaterial : public CMaterial
{
public:

	double Area;	//!< Sectional area of a bar element

public:
	
//!	����CBarMaterial�Ĳ��ϲ���
	virtual bool Read(ifstream& Input);

//!	���CBarMaterial�Ĳ��ϲ���
	virtual void Write(COutputter& output);
};

//! Q4���͵�Ԫ��Ӧ�Ĳ�����
class CQ4Material : public CMaterial
{
public:

	double Nu;	//!< Sectional Nu of a Q4 element

public:
	
//!	����CQ4Material�Ĳ��ϲ���
	virtual bool Read(ifstream& Input);

//!	���CQ4Material�Ĳ��ϲ���
	virtual void Write(COutputter& output);
};
