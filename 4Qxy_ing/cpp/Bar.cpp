/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Bar.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	构造函数，初始化各种东西
CBar::CBar()
{
	NEN_ = 2;	// Each element has 2 nodes
	nodes_ = new CNode*[NEN_];
    
    ND_ = 6;
	
    LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

//	析构函数
CBar::~CBar()
{
}

//  从ElementGroup对象中传入MaterialSets与NodeList，从文件中读入MSet与N1,N2,确定单元的节点以及所用的材料
bool CBar::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int MSet;	// Material property set number
	unsigned int N1, N2;	// Left node number and right node number

	Input >> N1 >> N2 >> MSet;
    //采用强制类型转换把指向CMaterial类型的指针转化为指向CBarMaterial类型的指针，并使用MSet获得指向该单元对应材料的指针
    ElementMaterial_ = dynamic_cast<CBarMaterial*>(MaterialSets) + MSet - 1;
	//从全局节点列表中提取该单元的节点并储存到nodes_数组里
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];

	return true;
}

//	输出Bar单元的属性
void CBar::Write(COutputter& output)
{
	//输出单元对应的两个节点的全局节点编号，单元
	output << setw(11) << nodes_[0]->NodeNumber
		   << setw(9) << nodes_[1]->NodeNumber << setw(12) << ElementMaterial_->nset << endl;
}

//	计算单元刚度阵并采用按列存储方法存到传入的Matrix数组里
void CBar::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());

//	算长度
	double DX[3];		//	dx = x2-x1, dy = y2-y1, dz = z2-z1
	for (unsigned int i = 0; i < 3; i++)
		DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];

	double DX2[6];	//  Quadratic polynomial (dx^2, dy^2, dz^2, dx*dy, dy*dz, dx*dz)
	DX2[0] = DX[0] * DX[0];
	DX2[1] = DX[1] * DX[1];
	DX2[2] = DX[2] * DX[2];
	DX2[3] = DX[0] * DX[1];
	DX2[4] = DX[1] * DX[2];
	DX2[5] = DX[0] * DX[2];

	double L2 = DX2[0] + DX2[1] + DX2[2];
	double L = sqrt(L2);

//	计算刚度阵，按列存进Matrix数组
//CBar类从CElement类继承了ElementMaterial_变量，这个变量的值在CBar::Read中确定
//这个变量是指向CMaterial类的指针，利用下面的语法可以把它改为指向CBarmaterial类的指针
	CBarMaterial* material_ = dynamic_cast<CBarMaterial*>(ElementMaterial_);	// Pointer to material of the element

	//用指向CBarmaterial类的指针读取各种材料参数
	double k = material_->E * material_->Area / L / L2;

	Matrix[0] = k*DX2[0];
	Matrix[1] = k*DX2[1];
	Matrix[2] = k*DX2[3];
	Matrix[3] = k*DX2[2];
	Matrix[4] = k*DX2[4];
	Matrix[5] = k*DX2[5];
	Matrix[6] = k*DX2[0];
	Matrix[7] = -k*DX2[5];
	Matrix[8] = -k*DX2[3];
	Matrix[9] = -k*DX2[0];
	Matrix[10] = k*DX2[1];
	Matrix[11] = k*DX2[3];
	Matrix[12] = -k*DX2[4];
	Matrix[13] = -k*DX2[1];
	Matrix[14] = -k*DX2[3];
	Matrix[15] = k*DX2[2];
	Matrix[16] = k*DX2[4];
	Matrix[17] = k*DX2[5];
	Matrix[18] = -k*DX2[2];
	Matrix[19] = -k*DX2[4];
	Matrix[20] = -k*DX2[5];
}

//	计算单元应力
void CBar::ElementStress(double* stress, double* Displacement)
{
	//把指向CMaterial类的指针转变成指向CBarMaterial类的指针
	CBarMaterial* material_ = dynamic_cast<CBarMaterial*>(ElementMaterial_);	

	double DX[3];	//	dx = x2-x1, dy = y2-y1, dz = z2-z1
	double L2 = 0;	//	Square of bar length (L^2)

	for (unsigned int i = 0; i < 3; i++)
	{
		DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];
		L2 = L2 + DX[i]*DX[i];
	}

	double S[6];
	for (unsigned int i = 0; i < 3; i++)
	{
		S[i] = -DX[i] * material_->E / L2;
		S[i+3] = -S[i];
	}
	//*代表对stress的解引用
	*stress = 0.0;
	for (unsigned int i = 0; i < 6; i++)
	{
		if (LocationMatrix_[i])
			*stress += S[i] * Displacement[LocationMatrix_[i]-1];
	}
}
