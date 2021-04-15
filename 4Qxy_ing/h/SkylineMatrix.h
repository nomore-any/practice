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

//! ���ڰ��д洢�ն�����࣬�ն���������������δ֪��ʹ��ģ����
template <class T_>
class CSkylineMatrix
{
    
private:
//! ����洢�����������ݵ�����
    T_* data_;
    
//! �ն����ά����Ҳ���Ƿ�����(number of equation)
    unsigned int NEQ_;

//! �������
    unsigned int MK_;

//! �洢�ն���������ڴ��С
    unsigned int NWK_;

//! �и�����
    unsigned int* ColumnHeights_;
    
//! �Խ�Ԫλ������
    unsigned int* DiagonalAddress_;
    
public:

//! ���캯�����ڸ�������(N=������)��û����ʱ���ò�ͬ�ĳ�ʼ������
    inline CSkylineMatrix();
    inline CSkylineMatrix(unsigned int N);
    
//! destructor
    inline ~CSkylineMatrix();

//! ����()��������Ա��þ������ʽ��������洢������(����K(1,2)�ڰ��д洢�������������ԭ��λ�ھ���(1,2)λ�õ���)
    inline T_& operator()(unsigned int i, unsigned int j);

#ifdef _DEBUG_
//! operator (i) where i numbering from 1
//! For the sake of efficiency, the index bounds are not checked
    inline T_& operator()(unsigned int i);
#endif
    
//! ���ն������洢�ռ�
    inline void Allocate();
    
//! ����ĳ����Ԫ�����LocationMatrix������õ�Ԫ���и�
    void CalculateColumnHeight(unsigned int* LocationMatrix, size_t ND);

//! �������� ( = max(ColumnHeights) + 1 ��
    void CalculateMaximumHalfBandwidth();

//! ����ն����жԽ�Ԫ�ĵ�ַ
    void CalculateDiagnoalAddress();

//! ����ĳ����Ԫ����������������Ѹõ�Ԫ�ĸն�����װ��ȫ�ָն���
    void Assembly(double* Matrix, unsigned int* LocationMatrix, size_t ND);

//! ����ָ���и������ָ��
    inline unsigned int* GetColumnHeights();

//! �����������
    inline unsigned int GetMaximumHalfBandwidth() const;

//! ����ָ��Խ�Ԫλ�������ָ��
    inline unsigned int* GetDiagonalAddress();

//! ���ؾ���ά��
    inline unsigned int dim() const;
    
//! ���ش洢�ն�������鳤��
    inline unsigned int size() const;

}; /* class definition */

//! ���캯��
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

//! ������������
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

//! ͨ������������ķ�ʽ���Ծ���ı�ʾ��ʽ(����K(1,2))���������ж�Ӧλ�õ���
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

//! Ϊ�ն����������洢�ռ�

template <class T_>
inline void CSkylineMatrix<T_>::Allocate()
{
    NWK_ = DiagonalAddress_[NEQ_] - DiagonalAddress_[0];
//��ע�⣬��CalculateDiagonalAddress()�����п��Կ���ʵ���϶Խ�Ԫλ���������NEQ_+1������������Ϊ��һ������NWK_�ṩ����
    
	//��������һ��NWK_���������ڴ洢�����ҳ�ʼ��
	data_ = new T_[NWK_];
    for (unsigned int i = 0; i < NWK_; i++)
        data_[i] = T_(0);
}

//! ����ָ���иߵ�ָ��
template <class T_>
inline unsigned int* CSkylineMatrix<T_>::GetColumnHeights()
{
    return ColumnHeights_;
}

//! �����������
template <class T_>
inline unsigned int CSkylineMatrix<T_>::GetMaximumHalfBandwidth() const
{
    return(MK_);
}

//! ����ָ��DiagonalAddress�����ָ��
template <class T_>
inline unsigned int* CSkylineMatrix<T_>::GetDiagonalAddress()
{
    return DiagonalAddress_;
}

//! ���ظն����ά��NEQ
template <class T_>
inline unsigned int CSkylineMatrix<T_>::dim() const
{
    return(NEQ_);
}

//! ���ش洢�ն�������鳤��
template <class T_>
inline unsigned int CSkylineMatrix<T_>::size() const
{
   return(NWK_);
}

//  ����һ����Ԫ���ж������Ԫ�Ƿ��ṩ������иߣ�����ǾͰ����и߸���ColumnHeights_[]
template <class T_>
void CSkylineMatrix<T_>::CalculateColumnHeight(unsigned int* LocationMatrix, size_t ND)
{    
//  Ѱ�Ҹõ�Ԫ���ɶȵ����ȫ�ֱ��
    unsigned int nfirstrow = INT_MAX;
    for (unsigned int i = 0; i < ND; i++)
        if (LocationMatrix[i] && LocationMatrix[i] < nfirstrow)
            nfirstrow = LocationMatrix[i];
    
//  ����õ�Ԫ��Ӧ������ȫ�����ɶȵ��иߣ���ColumnHeights_���Ƚ�
//  ĳ���иߵļ��㷽��Ϊ����Ԫ���ȫ�ֱ��-����
    for (unsigned int i = 0; i < ND; i++)
    {
        unsigned int column = LocationMatrix[i];
        if (!column)
            continue;
        
        unsigned int Height = column - nfirstrow;
        if (ColumnHeights_[column-1] < Height) ColumnHeights_[column-1] = Height;
    }
}

// ������� ( = max(ColumnHeights) + 1 ��
template <class T_>
void CSkylineMatrix<T_>::CalculateMaximumHalfBandwidth()
{
    MK_ = ColumnHeights_[0];
    
    for (unsigned int i=1; i<NEQ_; i++)
        if (MK_ < ColumnHeights_[i])
            MK_ = ColumnHeights_[i];
    
    MK_ = MK_ + 1;
}

//    ���뵥Ԫ�����������(�ն���(���д洢)��LocationMatrix��ND)���Ѹõ�Ԫ�ĸն�����װ��ȫ�ָն���
template <class T_>
void CSkylineMatrix<T_>::Assembly(double* Matrix, unsigned int* LocationMatrix, size_t ND)
{
//  �Ե�Ԫ�����ɶ�ѭ��(���ն�������б���)
    for (unsigned int j = 0; j < ND; j++)
    {
	//  ��ȡ��Ԫ�ն����j�ж�Ӧ��ȫ�ֱ��
        unsigned int Lj = LocationMatrix[j];    // Global equation number corresponding to jth DOF of the element
        if (!Lj) continue;
        
//      �����j�У���Ԫ�ն����жԽ�Ԫ�ڵ�Ԫ�ն������е�λ��
        unsigned int DiagjElement = (j+1)*j/2;
        
		//�ڸն���ĵ�j�ж���ѭ��
        for (unsigned int i = 0; i <= j; i++)
        {   
			//��ȡ��i�ж�Ӧ��ȫ�ֱ��
            unsigned int Li = LocationMatrix[i];    // Global equation number corresponding to ith DOF of the element
            
            if (!Li) continue;
            //ʹ����SkylineMatrix�����ص�()��������ѵ�Ԫ�ն�������ݼӽ�ȫ�ָն���
            (*this)(Li,Lj) += Matrix[DiagjElement + j - i];
        }
    }
    
    return;
}

//    ����ն���Խ�Ԫ�������е�λ��
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

