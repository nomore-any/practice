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

//!	���һ������
template <class type> void clear( type* a, unsigned int N );

//!	CDomain�ඨ�������
/*!	Only a single instance of Domain class can be created */
class CDomain
{
private:

//!	����һ��ָ��CDomain���ָ�룬���ָ���Ǿ�̬�ģ�������������ԣ��������κ�һ�����������
	static CDomain* _instance;

//!	����һ��ifstream��Ķ����ṩ�����ļ��Ľӿڣ���������������붼Ҫ����Input
	ifstream Input;

//!	����
	char Title[256]; 

//!	�������ģʽ
/*!		0 : Data check only;
		1 : Execution */
	unsigned int MODEX;

//!	�ڵ�����
	unsigned int NUMNP;

//!	�ڵ��б�(�ڵ����͵�����)
	CNode* NodeList;

//!	��Ԫ������
//! һ����Ԫ��ֻ��һ�����͵ĵ�Ԫ
	unsigned int NUMEG;

//! ��Ԫ���б�(��Ԫ�����͵�����)
    CElementGroup* EleGrpList;
    
//!	�غɹ�����
	unsigned int NLCASE;

//!	�غɹ����б�(�غɹ����������͵����飩
	CLoadCaseData* LoadCases;

//!	�غ����б�(ÿ���غɹ����ж����غ���)
	unsigned int* NLOAD;

//!	�ܷ�����(��δԼ�����ɶ�����
	unsigned int NEQ;

//!	����һ������ȫ�ָն��������
/*! A one-dimensional array storing only the elements below the	skyline of the 
    global stiffness matrix. */
    CSkylineMatrix<double>* StiffnessMatrix;

//!	����һ������ȫ���غ�����������
	double* Force;

private:

//!	����
	CDomain();

//!	����
	~CDomain();

public:

//!	����һ��ָ��CDomain�����ָ�룬�������������CDomain�����������
	static CDomain* GetInstance();

//!	�����ļ���
	bool ReadData(string FileName, string OutFile);

//!	��ȡ�ڵ�����
	bool ReadNodalPoints();

//!	��ȡ�غ�����
	bool ReadLoadCases();

//!	��ȡ��Ԫ����
	bool ReadElements();

//!	���㷽����NEQ��˳���ÿ����Ԫ����LocationMatrix
	void CalculateEquationNumber();

//!	����ȫ�ָն����и�
	void CalculateColumnHeights();

//! ���ն��������������ڴ棬������һЩ�ն������
/*!	Allocate Force, ColumnHeights, DiagonalAddress and StiffnessMatrix and 
    calculate the column heights and address of diagonal elements */
	void AllocateMatrices();

//!	��װȫ�ָն���
	void AssembleStiffnessMatrix();

//!	��LoadCase�غɹ�����װ��ȫ����������
	bool AssembleForce(unsigned int LoadCase); 

	//�����Ƿ��ظ�������

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
