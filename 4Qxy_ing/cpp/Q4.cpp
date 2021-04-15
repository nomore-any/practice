#include "Q4.h"
#include "math.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	���캯������ʼ�����ֶ���
CQ4::CQ4()
{
	NEN_ = 4;	// Each element has 4 nodes
	nodes_ = new CNode*[NEN_];
    
    ND_ = 8;
	
    LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

//	��������
CQ4::~CQ4()
{
}

//  ��ElementGroup�����д���MaterialSets��NodeList�����ļ��ж���MSet��N1,N2,N3,N4,ȷ����Ԫ�Ľڵ��Լ����õĲ���
bool CQ4::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int MSet;	// Material property set number
	unsigned int N1, N2, N3, N4;	// Left node number and right node number

	Input >> N1 >> N2 >> N3 >> N4 >> MSet;
    //����ǿ������ת����ָ��CMaterial���͵�ָ��ת��Ϊָ��CQ4Material���͵�ָ�룬��ʹ��MSet���ָ��õ�Ԫ��Ӧ���ϵ�ָ��
    ElementMaterial_ = dynamic_cast<CQ4Material*>(MaterialSets) + MSet - 1;
	//��ȫ�ֽڵ��б�����ȡ�õ�Ԫ�Ľڵ㲢���浽nodes_������
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
    nodes_[2] = &NodeList[N3 - 1];
	nodes_[3] = &NodeList[N4 - 1];
	return true;
}

//	���Q4CQ4��Ԫ������
void CQ4::Write(COutputter& output)
{
	//�����Ԫ��Ӧ���ĸ��ڵ��ȫ�ֽڵ��ţ���Ԫ
	output << setw(11) << nodes_[0]->NodeNumber
		   << setw(9) << nodes_[1]->NodeNumber 
           << setw(9) << nodes_[2]->NodeNumber
		   << setw(9) << nodes_[3]->NodeNumber << setw(12) << ElementMaterial_->nset << endl;
}

//	���㵥Ԫ�ն��󲢲��ð��д洢�����浽�����matrix������
void CQ4::ElementStiffness(double* matrix)
{
	clear(matrix, SizeOfStiffnessMatrix());

//	�����꣬������Ӧ��
	double x1,x2,x3,x4,y1,y2,y3,y4;
 	x1= nodes_[0]->XYZ[0];
 	x2= nodes_[1]->XYZ[0];
	x3= nodes_[2]->XYZ[0];
	x4= nodes_[3]->XYZ[0];
	y1= nodes_[0]->XYZ[1];
	y2= nodes_[1]->XYZ[2];
	y3= nodes_[2]->XYZ[3];
	y4= nodes_[3]->XYZ[4];

	double x21 = x2-x1;
	double x43 = x4-x3;
	double y21 = y2-y1;
	double y43 = y4-y3;
	double x41 = x4-x1;
	double x32 = x3-x2;
	double y41 = y4-y1;
	double y32 = y3-y2;
	double f[2],g[2];
	f[0] = sqrt(3)/3;
	f[1] = -sqrt(3)/3;
	g[0] = sqrt(3)/3;
	g[1] = -sqrt(3)/3;

//	����ն��󣬰��д��matrix����

//CQ4���CElement��̳���ElementMaterial_���������������ֵ��CQ4::Read��ȷ��
//���������ָ��CMaterial���ָ�룬����������﷨���԰�����Ϊָ��CQ4material���ָ��
	CQ4Material* material_ = dynamic_cast<CQ4Material*>(ElementMaterial_);	// Pointer to material of the element

	//��ָ��CQ4material���ָ���ȡ���ֲ��ϲ���
	double k = material_->E / (1-(material_-> Nu)*(material_-> Nu));
	double u = material_-> Nu;


	for(int i = 0; i < 2; i++)
	{	for(int j = 0; j < 2; j++) 
		{	double J1[2][2],J2[2][2],J3[2][2],J4[2][2];
			J1[i][j] = 0.25*(x21*(1-f[i])-x43*(1+f[i]));
			J2[i][j] = 0.25*(y21*(1-f[i])-y43*(1+f[i]));
			J3[i][j] = 0.25*(x41*(1-g[j])-x32*(1+g[j]));
			J4[i][j] = 0.25*(y41*(1-g[j])-y32*(1+g[j]));
			
			double Nx1[i][j],Nx2[i][j],Nx3[i][j],Nx4[i][j],Ny1[i][j],Ny2[i][j],Ny3[i][j],Ny4[i][j];
			Nx1[i][j] = 0.25*(J4[i][j]*((f[i])-1)-J2[i][j]*(g[j]-1));
			Nx2[i][j] = 0.25*(J4[i][j]*(1-(f[i]))-J2[i][j]*(-g[j]-1));
			Nx3[i][j] = 0.25*(J4[i][j]*(1+(f[i]))-J2[i][j]*(g[j]+1));
			Nx4[i][j] = 0.25*(J4[i][j]*(-(f[i])-1)-J2[i][j]*(-g[j]+1));
			Nx1[i][j] = 0.25*(-J3[i][j]*((f[i])-1)+J1[i][j]*(g[j]-1));
			Nx2[i][j] = 0.25*(-J3[i][j]*(1-(f[i]))+J1[i][j]*(-g[j]-1));
			Nx3[i][j] = 0.25*(-J3[i][j]*(1+(f[i]))+J1[i][j]*(g[j]+1));
			Nx4[i][j] = 0.25*(-J3[i][j]*(-(f[i])-1)+J1[i][j]*(-g[j]+1));
			
			double detJ[2][2];
			detJ[i][j] = J1[i][j]*J4[i][j]-J2[i][j]*J3[i][j];

			double matrix[36];
			matrix[0]=k*(Nx1[i][j]*Nx1[i][j]+(0.5-0.5*u)*Ny1[i][j]*Ny1[i][j])/detJ[i][j];
			matrix[1]=k*(Ny1[i][j]*Ny1[i][j]+(0.5-0.5*u)*Nx1[i][j]*Nx1[i][j])/detJ[i][j];
			matrix[2]=k*(Ny1[i][j]*Nx1[i][j]*(0.5+0.5*u))/detJ[i][j];
			matrix[3]=k*(Nx2[i][j]*Nx2[i][j]+(0.5-0.5*u)*Ny2[i][j]*Ny2[i][j])/detJ[i][j];
			matrix[4]=k*(u*Ny1[i][j]*Nx2[i][j]+(0.5-0.5*u)*Ny2[i][j]*Nx1[i][j])/detJ[i][j];
			matrix[5]=k*(Nx2[i][j]*Nx1[i][j]+(0.5-0.5*u)*Ny2[i][j]*Ny1[i][j])/detJ[i][j];
			matrix[6]=k*(Ny2[i][j]*Ny2[i][j]+(0.5-0.5*u)*Nx2[i][j]*Nx2[i][j])/detJ[i][j];
			matrix[7]=k*(Ny2[i][j]*Nx2[i][j]*(0.5+0.5*u))/detJ[i][j];
			matrix[8]=k*(Ny1[i][j]*Ny2[i][j]+(0.5-0.5*u)*Nx1[i][j]*Nx2[i][j])/detJ[i][j];
			matrix[9]=k*(u*Ny2[i][j]*Nx1[i][j]+(0.5-0.5*u)*Nx2[i][j]*Ny1[i][j])/detJ[i][j];
			matrix[10]=k*(Nx3[i][j]*Nx3[i][j]+(0.5-0.5*u)*Ny3[i][j]*Ny3[i][j])/detJ[i][j];
			matrix[11]=k*(u*Ny2[i][j]*Nx3[i][j]+(0.5-0.5*u)*Nx2[i][j]*Ny3[i][j])/detJ[i][j];
			matrix[12]=k*(Nx2[i][j]*Nx3[i][j]+(0.5-0.5*u)*Ny2[i][j]*Ny3[i][j])/detJ[i][j];
			matrix[13]=k*(u*Ny1[i][j]*Nx3[i][j]+(0.5-0.5*u)*Nx1[i][j]*Ny3[i][j])/detJ[i][j];
			matrix[14]=k*(Nx1[i][j]*Nx3[i][j]+(0.5-0.5*u)*Ny1[i][j]*Ny3[i][j])/detJ[i][j];
			matrix[15]=k*(Ny3[i][j]*Ny3[i][j]+(0.5-0.5*u)*Nx3[i][j]*Nx3[i][j])/detJ[i][j];
			matrix[16]=k*((0.5+0.5*u)*Nx3[i][j]*Ny3[i][j])/detJ[i][j];
			matrix[17]=k*(Ny2[i][j]*Ny3[i][j]+(0.5-0.5*u)*Nx2[i][j]*Nx3[i][j])/detJ[i][j];
			matrix[18]=k*(u*Nx2[i][j]*Ny3[i][j]+(0.5-0.5*u)*Ny2[i][j]*Nx3[i][j])/detJ[i][j];
			matrix[19]=k*(Ny1[i][j]*Ny3[i][j]+(0.5-0.5*u)*Nx1[i][j]*Nx3[i][j])/detJ[i][j];
			matrix[20]=k*(u*Nx1[i][j]*Ny3[i][j]+(0.5-0.5*u)*Ny1[i][j]*Nx3[i][j])/detJ[i][j];
			matrix[21]=k*(Nx4[i][j]*Nx4[i][j]+(0.5-0.5*u)*Ny4[i][j]*Ny4[i][j])/detJ[i][j];
			matrix[22]=k*(u*Ny3[i][j]*Nx4[i][j]+(0.5-0.5*u)*Nx3[i][j]*Ny4[i][j])/detJ[i][j];
			matrix[23]=k*(Nx3[i][j]*Nx4[i][j]+(0.5-0.5*u)*Ny3[i][j]*Ny4[i][j])/detJ[i][j];
			matrix[24]=k*(u*Ny2[i][j]*Nx4[i][j]+(0.5-0.5*u)*Nx2[i][j]*Ny4[i][j])/detJ[i][j];
			matrix[25]=k*(Nx2[i][j]*Nx4[i][j]+(0.5-0.5*u)*Ny2[i][j]*Ny4[i][j])/detJ[i][j];
			matrix[26]=k*(u*Ny1[i][j]*Nx4[i][j]+(0.5-0.5*u)*Nx1[i][j]*Ny4[i][j])/detJ[i][j];
			matrix[27]=k*(Nx1[i][j]*Nx4[i][j]+(0.5-0.5*u)*Ny1[i][j]*Ny4[i][j])/detJ[i][j];
			matrix[28]=k*(Ny4[i][j]*Ny4[i][j]+(0.5-0.5*u)*Nx4[i][j]*Nx4[i][j])/detJ[i][j];
			matrix[29]=k*(u*Nx4[i][j]*Ny4[i][j]+(0.5-0.5*u)*Ny4[i][j]*Nx4[i][j])/detJ[i][j];
			matrix[30]=k*(Ny3[i][j]*Ny4[i][j]+(0.5-0.5*u)*Nx3[i][j]*Nx4[i][j])/detJ[i][j];
			matrix[31]=k*(u*Nx3[i][j]*Ny4[i][j]+(0.5-0.5*u)*Ny3[i][j]*Nx4[i][j])/detJ[i][j];
			matrix[32]=k*(Ny2[i][j]*Ny4[i][j]+(0.5-0.5*u)*Nx2[i][j]*Nx4[i][j])/detJ[i][j];
			matrix[33]=k*(u*Nx2[i][j]*Ny4[i][j]+(0.5-0.5*u)*Ny2[i][j]*Nx4[i][j])/detJ[i][j];
			matrix[34]=k*(Ny1[i][j]*Ny4[i][j]+(0.5-0.5*u)*Nx1[i][j]*Nx4[i][j])/detJ[i][j];
			matrix[35]=k*(u*Nx1[i][j]*Ny4[i][j]+(0.5-0.5*u)*Ny1[i][j]*Nx4[i][j])/detJ[i][j];
			for(int m = 0; m < 36; m++)
			{
				double Matrix[36];
				Matrix[m]=matrix[m]+Matrix[m];
			}
		}
	}
}

//	���㵥ԪӦ��
/*void CQ4::ElementStress(double* stress, double* Displacement)
{
	//��ָ��CMaterial���ָ��ת���ָ��CQ4Material���ָ��
	CQ4Material* material_ = dynamic_cast<CQ4Material*>(ElementMaterial_);	

	double DX[3];	//	dx = x2-x1, dy = y2-y1, dz = z2-z1
	double L2 = 0;	//	Square of Q4CQ4 length (L^2)

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
	//*�����stress�Ľ�����
	*stress = 0.0;
	for (unsigned int i = 0; i < 6; i++)
	{
		if (Locationmatrix_[i])
			*stress += S[i] * Displacement[Locationmatrix_[i]-1];
	}
}
*/

//����ն���
	
