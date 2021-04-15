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

#include "Node.h"
#include "ElementGroup.h"
#include "Outputter.h"
#include "Solver.h"
#include "LoadCaseData.h"
#include "SkylineMatrix.h"

using namespace std;

//!	清除一个数组
template <class type> void clear( type* a, unsigned int N );

//!	CDomain类定义求解域
/*!	Only a single instance of Domain class can be created */
class CDomain
{
private:

//!	创建一个指向CDomain类的指针，这个指针是静态的，所以是类的属性，而不是任何一个对象的属性
	static CDomain* _instance;

//!	创建一个ifstream类的对象，提供输入文件的接口，所有其他类的输入都要借助Input
	ifstream Input;

//!	标题
	char Title[256]; 

//!	定义求解模式
/*!		0 : Data check only;
		1 : Execution */
	unsigned int MODEX;

//!	节点总数
	unsigned int NUMNP;

//!	节点列表(节点类型的数组)
	CNode* NodeList;

//!	单元组总数
//! 一个单元组只有一种类型的单元
	unsigned int NUMEG;

//! 单元组列表(单元组类型的数组)
    CElementGroup* EleGrpList;
    
//!	载荷工况数
	unsigned int NLCASE;

//!	载荷工况列表(载荷工况数据类型的数组）
	CLoadCaseData* LoadCases;

//!	载荷数列表(每种载荷工况有多少载荷数)
	unsigned int* NLOAD;

//!	总方程数(总未约束自由度数）
	unsigned int NEQ;

//!	定义一个储存全局刚度阵的数组
/*! A one-dimensional array storing only the elements below the	skyline of the 
    global stiffness matrix. */
    CSkylineMatrix<double>* StiffnessMatrix;

//!	定义一个储存全局载荷向量的数组
	double* Force;

private:

//!	构造
	CDomain();

//!	析构
	~CDomain();

public:

//!	返回一个指向CDomain对象的指针，便于其他类调用CDomain对象里的数据
	static CDomain* GetInstance();

//!	输入文件名
	bool ReadData(string FileName, string OutFile);

//!	读取节点数据
	bool ReadNodalPoints();

//!	读取载荷数据
	bool ReadLoadCases();

//!	读取单元数据
	bool ReadElements();

//!	计算方程数NEQ，顺便给每个单元生成LocationMatrix
	void CalculateEquationNumber();

//!	计算全局刚度阵列高
	void CalculateColumnHeights();

//! 给刚度阵、力向量分配内存，并计算一些刚度阵变量
/*!	Allocate Force, ColumnHeights, DiagonalAddress and StiffnessMatrix and 
    calculate the column heights and address of diagonal elements */
	void AllocateMatrices();

//!	组装全局刚度阵
	void AssembleStiffnessMatrix();

//!	把LoadCase载荷工况组装到全局力向量里
	bool AssembleForce(unsigned int LoadCase); 

	//下面是返回各种属性

//!	Return solution mode
	inline unsigned int GetMODEX() { return MODEX; }

//!	Return the title of problem
	inline string GetTitle() { return Title; }

//!	Return the total number of equations
	inline unsigned int GetNEQ() { return NEQ; }

//!	Return the total number of nodal points
	inline unsigned int GetNUMNP() { return NUMNP; }

//!	Return the node list
	inline CNode* GetNodeList() { return NodeList; }

//!	Return total number of element groups
	inline unsigned int GetNUMEG() { return NUMEG; }

//! Return element group list
    inline CElementGroup* GetEleGrpList() { return EleGrpList; }

//!	Return pointer to the global nodal force vector
	inline double* GetForce() { return Force; }

//!	Return pointer to the global nodal displacement vector
	inline double* GetDisplacement() { return Force; }

//!	Return the total number of load cases
	inline unsigned int GetNLCASE() { return NLCASE; }

//!	Return the number of concentrated loads applied in each load case
	inline unsigned int* GetNLOAD() { return NLOAD; }

//!	Return the list of load cases
	inline CLoadCaseData* GetLoadCases() { return LoadCases; }

//!	Return pointer to the banded stiffness matrix
	inline CSkylineMatrix<double>* GetStiffnessMatrix() { return StiffnessMatrix; }

};
