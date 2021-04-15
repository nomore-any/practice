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
#include "Material.h"

using namespace std;

template <class type> void clear( type* a, unsigned int N );	// Clear an array

//������һ����Ԫ�ĸ��࣬���������������е�Ԫ�����е�����
class CElement
{
protected:

//!	����ÿ����Ԫ�Ľڵ���
	unsigned int NEN_;

//!	����һ����ָ��ڵ��ָ�빹�ɵ�����
	CNode** nodes_;

//!	����һ��ָ��CMaterial���ָ��
	CMaterial* ElementMaterial_;	//!< Pointer to an element of MaterialSetList[][]
    
//! ����������(��һ����Ԫ��˵��һ�����飬��ӳ�ֲ���ź�ȫ�ֱ�Ź�ϵ)
    unsigned int* LocationMatrix_;

//! LocationMatrix�ĳ��ȣ���ʵҲ����һ����Ԫ�������ɶȸ���
    unsigned int ND_;

public:

//!	���캯�����Ը��ֶ������г�ʼ��
	CElement() : NEN_(0), nodes_(nullptr), ElementMaterial_(nullptr) {}

//! ���������������ù�
    virtual ~CElement() {
        if (!nodes_)
            delete [] nodes_;
        
        if (!ElementMaterial_)
            delete [] ElementMaterial_;
        
        if (!LocationMatrix_)
            delete [] LocationMatrix_;
    }

//!	�麯������ȡ��Ԫ���ݣ�����������������Ϊ����ElementGroup�����涨���Read�����ܹ��ڶ�ȡ��ͬʱ˳���ж������ʽ�Ƿ���ȷ
	virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList) = 0;

//!	�����Ԫ������
	virtual void Write(COutputter& output) = 0;

//! �����������ʱ��CNode���bcode�������Ѿ������˽ڵ��Ӧ���ɶȵ�ȫ�ֱ�ţ����Դ�����ȡ���õ�Ԫ��LocationMatrix
    virtual void GenerateLocationMatrix()
    {
        unsigned int i = 0;
        for (unsigned int N = 0; N < NEN_; N++)
            for (unsigned int D = 0; D < CNode::NDF; D++)
                LocationMatrix_[i++] = nodes_[N]->bcode[D];
    }

//! ���㵥Ԫ�ն�����ð��д洢��������Ҫ�������С
    virtual unsigned int SizeOfStiffnessMatrix()
    {
        unsigned int size = 0;
        for (int i=1; i<= ND_; i++)
            size += i;
        
        return size;
    }

//!	���㵥Ԫ�ն���(���д洢��������ʽ)���洢�������stiffness������
	virtual void ElementStiffness(double* stiffness) = 0; 

//!	���ݴ����ȫ�ֽڵ�λ�ƣ����㵥Ԫ��Ӧ��
	virtual void ElementStress(double* stress, double* Displacement) = 0;

//! ����һ����Ԫ�Ľڵ���
    inline unsigned int GetNEN() { return NEN_; }
    
//!	����nodes_
	inline CNode** GetNodes() { return nodes_; }

//!	���ص�Ԫ���ϵ�ָ�룬����һ��ָ��CMaterial���͵�ָ��
	inline CMaterial* GetElementMaterial() { return ElementMaterial_; }
    
    //! ���ص�Ԫ�ķ�������
    inline unsigned int* GetLocationMatrix() { return LocationMatrix_; }
    
    //! ���ص�Ԫ�����ɶ���
    inline unsigned int GetND() { return ND_; }
};
