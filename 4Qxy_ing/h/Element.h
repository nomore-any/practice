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

//定义了一个单元的父类，这个父类包含了所有单元都具有的属性
class CElement
{
protected:

//!	定义每个单元的节点数
	unsigned int NEN_;

//!	定义一个由指向节点的指针构成的数组
	CNode** nodes_;

//!	定义一个指向CMaterial类的指针
	CMaterial* ElementMaterial_;	//!< Pointer to an element of MaterialSetList[][]
    
//! 定义分配矩阵(对一个单元来说是一个数组，反映局部编号和全局编号关系)
    unsigned int* LocationMatrix_;

//! LocationMatrix的长度，其实也就是一个单元的总自由度个数
    unsigned int ND_;

public:

//!	构造函数，对各种东西进行初始化
	CElement() : NEN_(0), nodes_(nullptr), ElementMaterial_(nullptr) {}

//! 虚析构函数，不用管
    virtual ~CElement() {
        if (!nodes_)
            delete [] nodes_;
        
        if (!ElementMaterial_)
            delete [] ElementMaterial_;
        
        if (!LocationMatrix_)
            delete [] LocationMatrix_;
    }

//!	虚函数，读取单元数据，布尔运算符的输出是为了在ElementGroup类里面定义的Read函数能够在读取的同时顺便判断输入格式是否正确
	virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList) = 0;

//!	输出单元的数据
	virtual void Write(COutputter& output) = 0;

//! 调用这个函数时，CNode类的bcode数组中已经包含了节点对应自由度的全局编号，可以从中提取出该单元的LocationMatrix
    virtual void GenerateLocationMatrix()
    {
        unsigned int i = 0;
        for (unsigned int N = 0; N < NEN_; N++)
            for (unsigned int D = 0; D < CNode::NDF; D++)
                LocationMatrix_[i++] = nodes_[N]->bcode[D];
    }

//! 计算单元刚度阵采用按列存储方法所需要的数组大小
    virtual unsigned int SizeOfStiffnessMatrix()
    {
        unsigned int size = 0;
        for (int i=1; i<= ND_; i++)
            size += i;
        
        return size;
    }

//!	计算单元刚度阵(按列存储，矩阵形式)，存储到传入的stiffness数组中
	virtual void ElementStiffness(double* stiffness) = 0; 

//!	根据传入的全局节点位移，计算单元的应力
	virtual void ElementStress(double* stress, double* Displacement) = 0;

//! 返回一个单元的节点数
    inline unsigned int GetNEN() { return NEN_; }
    
//!	返回nodes_
	inline CNode** GetNodes() { return nodes_; }

//!	返回单元材料的指针，这是一个指向CMaterial类型的指针
	inline CMaterial* GetElementMaterial() { return ElementMaterial_; }
    
    //! 返回单元的分配数组
    inline unsigned int* GetLocationMatrix() { return LocationMatrix_; }
    
    //! 返回单元的自由度数
    inline unsigned int GetND() { return ND_; }
};
