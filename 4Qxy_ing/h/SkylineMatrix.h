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

#include <string>
#include <climits>

#ifdef _DEBUG_
#include "Outputter.h"
#endif

//! 用于按列存储刚度阵的类，刚度阵数据类型事先未知，使用模板类
template <class T_>
class CSkylineMatrix
{
    
private:
//! 定义存储任意类型数据的数组
    T_* data_;
    
//! 刚度阵的维数，也就是方程数(number of equation)
    unsigned int NEQ_;

//! 最大半带宽
    unsigned int MK_;

//! 存储刚度阵所需的内存大小
    unsigned int NWK_;

//! 列高数组
    unsigned int* ColumnHeights_;
    
//! 对角元位置数组
    unsigned int* DiagonalAddress_;
    
public:

//! 构造函数，在给定输入(N=方程数)和没给定时采用不同的初始化方法
    inline CSkylineMatrix();
    inline CSkylineMatrix(unsigned int N);
    
//! destructor
    inline ~CSkylineMatrix();

//! 重载()运算符，以便用矩阵的形式调用数组存储的内容(比如K(1,2)在按列存储的数组里面调用原本位于矩阵(1,2)位置的数)
    inline T_& operator()(unsigned int i, unsigned int j);

#ifdef _DEBUG_
//! operator (i) where i numbering from 1
//! For the sake of efficiency, the index bounds are not checked
    inline T_& operator()(unsigned int i);
#endif
    
//! 给刚度阵分配存储空间
    inline void Allocate();
    
//! 传入某个单元对象的LocationMatrix，计算该单元的列高
    void CalculateColumnHeight(unsigned int* LocationMatrix, size_t ND);

//! 计算半带宽 ( = max(ColumnHeights) + 1 ）
    void CalculateMaximumHalfBandwidth();

//! 计算刚度阵中对角元的地址
    void CalculateDiagnoalAddress();

//! 传入某个单元对象的三个参数，把该单元的刚度阵组装进全局刚度阵
    void Assembly(double* Matrix, unsigned int* LocationMatrix, size_t ND);

//! 返回指向列高数组的指针
    inline unsigned int* GetColumnHeights();

//! 返回最大半带宽
    inline unsigned int GetMaximumHalfBandwidth() const;

//! 返回指向对角元位置数组的指针
    inline unsigned int* GetDiagonalAddress();

//! 返回矩阵维数
    inline unsigned int dim() const;
    
//! 返回存储刚度阵的数组长度
    inline unsigned int size() const;

}; /* class definition */

//! 构造函数
template <class T_>
inline CSkylineMatrix<T_>::CSkylineMatrix()
{
    NEQ_ = 0;
    MK_  = 0;
    NWK_ = 0;
    
    data_ = nullptr;
    ColumnHeights_ = nullptr;
    DiagonalAddress_ = nullptr;
}

template <class T_>
inline CSkylineMatrix<T_>::CSkylineMatrix(unsigned int N)
{
    NEQ_ = N;
    MK_  = 0;
    NWK_ = 0;

    data_ = nullptr;
    
    ColumnHeights_ = new unsigned int [NEQ_];
    for (unsigned int i = 0; i < NEQ_; i++)
        ColumnHeights_[i] = 0;

    DiagonalAddress_ = new unsigned int [NEQ_ + 1];
    for (unsigned int i = 0; i < NEQ_ + 1; i++)
        DiagonalAddress_[i] = 0;
}

//! 析构函数不管
template <class T_>
inline CSkylineMatrix<T_>::~CSkylineMatrix<T_>()
{
    if (ColumnHeights_)
        delete[] ColumnHeights_;
    
    if (DiagonalAddress_)
        delete[] DiagonalAddress_;
    
    if (data_)
        delete[] data_;
}

//! 通过重载运算符的方式，以矩阵的表示形式(比如K(1,2))调用数组中对应位置的数
template <class T_>
inline T_& CSkylineMatrix<T_>::operator()(unsigned int i, unsigned int j)
{
    if (j >= i)
        return data_[DiagonalAddress_[j - 1] + (j - i) - 1];
    else
        return data_[DiagonalAddress_[i - 1] + (i - j) - 1];
}

#ifdef _DEBUG_
//! operator function (i) where i numbering from 1
template <class T_>
inline T_& CSkylineMatrix<T_>::operator()(unsigned int i)
{
    return data_[i];
}
#endif

//! 为刚度阵数组分配存储空间

template <class T_>
inline void CSkylineMatrix<T_>::Allocate()
{
    NWK_ = DiagonalAddress_[NEQ_] - DiagonalAddress_[0];
//！注意，从CalculateDiagonalAddress()函数中可以看出实际上对角元位置数组存了NEQ_+1个数，就是在为这一步计算NWK_提供便利
    
	//堆区开辟一个NWK_的数组用于存储，并且初始化
	data_ = new T_[NWK_];
    for (unsigned int i = 0; i < NWK_; i++)
        data_[i] = T_(0);
}

//! 返回指向列高的指针
template <class T_>
inline unsigned int* CSkylineMatrix<T_>::GetColumnHeights()
{
    return ColumnHeights_;
}

//! 返回最大半带宽
template <class T_>
inline unsigned int CSkylineMatrix<T_>::GetMaximumHalfBandwidth() const
{
    return(MK_);
}

//! 返回指向DiagonalAddress数组的指针
template <class T_>
inline unsigned int* CSkylineMatrix<T_>::GetDiagonalAddress()
{
    return DiagonalAddress_;
}

//! 返回刚度阵的维数NEQ
template <class T_>
inline unsigned int CSkylineMatrix<T_>::dim() const
{
    return(NEQ_);
}

//! 返回存储刚度阵的数组长度
template <class T_>
inline unsigned int CSkylineMatrix<T_>::size() const
{
   return(NWK_);
}

//  给定一个单元，判断这个单元是否提供更大的列高，如果是就把新列高赋给ColumnHeights_[]
template <class T_>
void CSkylineMatrix<T_>::CalculateColumnHeight(unsigned int* LocationMatrix, size_t ND)
{    
//  寻找该单元自由度的最大全局编号
    unsigned int nfirstrow = INT_MAX;
    for (unsigned int i = 0; i < ND; i++)
        if (LocationMatrix[i] && LocationMatrix[i] < nfirstrow)
            nfirstrow = LocationMatrix[i];
    
//  计算该单元对应的所有全局自由度的列高，与ColumnHeights_作比较
//  某列列高的计算方法为：单元最大全局编号-列数
    for (unsigned int i = 0; i < ND; i++)
    {
        unsigned int column = LocationMatrix[i];
        if (!column)
            continue;
        
        unsigned int Height = column - nfirstrow;
        if (ColumnHeights_[column-1] < Height) ColumnHeights_[column-1] = Height;
    }
}

// 最大半带宽 ( = max(ColumnHeights) + 1 ）
template <class T_>
void CSkylineMatrix<T_>::CalculateMaximumHalfBandwidth()
{
    MK_ = ColumnHeights_[0];
    
    for (unsigned int i=1; i<NEQ_; i++)
        if (MK_ < ColumnHeights_[i])
            MK_ = ColumnHeights_[i];
    
    MK_ = MK_ + 1;
}

//    传入单元类的三个变量(刚度阵(按列存储)，LocationMatrix，ND)，把该单元的刚度阵组装到全局刚度阵
template <class T_>
void CSkylineMatrix<T_>::Assembly(double* Matrix, unsigned int* LocationMatrix, size_t ND)
{
//  对单元的自由度循环(即刚度阵里对列遍历)
    for (unsigned int j = 0; j < ND; j++)
    {
	//  获取单元刚度阵第j列对应的全局编号
        unsigned int Lj = LocationMatrix[j];    // Global equation number corresponding to jth DOF of the element
        if (!Lj) continue;
        
//      计算第j列，单元刚度阵中对角元在单元刚度数组中的位置
        unsigned int DiagjElement = (j+1)*j/2;
        
		//在刚度阵的第j列对行循环
        for (unsigned int i = 0; i <= j; i++)
        {   
			//获取第i行对应的全局编号
            unsigned int Li = LocationMatrix[i];    // Global equation number corresponding to ith DOF of the element
            
            if (!Li) continue;
            //使用了SkylineMatrix类重载的()运算符，把单元刚度阵的数据加进全局刚度阵
            (*this)(Li,Lj) += Matrix[DiagjElement + j - i];
        }
    }
    
    return;
}

//    计算刚度阵对角元在数组中的位置
template <class T_>
void CSkylineMatrix<T_>::CalculateDiagnoalAddress()
{
    //    Calculate the address of diagonal elements
    //    M(0) = 1;  M(i+1) = M(i) + H(i) + 1 (i = 0:NEQ)
    DiagonalAddress_[0] = 1;
    for (unsigned int col = 1; col <= NEQ_; col++)
        DiagonalAddress_[col] = DiagonalAddress_[col - 1] + ColumnHeights_[col-1] + 1;
    
#ifdef _DEBUG_
    COutputter* Output = COutputter::GetInstance();
    Output->PrintDiagonalAddress();
#endif
    
}

