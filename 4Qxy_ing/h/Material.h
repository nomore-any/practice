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
//定义一个材料的父类，这个父类中定义了所有材料都具有的属性
//!	Material base class which only define one data member
/*!	All type of material classes should be derived from this base class */
class CMaterial
{
public:

	unsigned int nset;	//！定义材料的编号，材料编号1,2,3...是随ElementGroup变化而刷新的
	
	double E;  //!定义材料的杨氏模量

public:

//! 虚析构函数
    virtual ~CMaterial() {};

//!	虚函数，用于读取材料参数
	virtual bool Read(ifstream& Input) = 0;

//!	虚函数，用于输出材料参数
    virtual void Write(COutputter& output) = 0;

};

//!	Bar类型单元对应的材料类
class CBarMaterial : public CMaterial
{
public:

	double Area;	//!< Sectional area of a bar element

public:
	
//!	读入CBarMaterial的材料参数
	virtual bool Read(ifstream& Input);

//!	输出CBarMaterial的材料参数
	virtual void Write(COutputter& output);
};

//! Q4类型单元对应的材料类
class CQ4Material : public CMaterial
{
public:

	double Nu;	//!< Sectional Nu of a Q4 element

public:
	
//!	读入CQ4Material的材料参数
	virtual bool Read(ifstream& Input);

//!	输出CQ4Material的材料参数
	virtual void Write(COutputter& output);
};
