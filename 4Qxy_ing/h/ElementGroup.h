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

#include <fstream>

#include "Element.h"
#include "Bar.h"
#include "Q4.h"
#include "Material.h"
#include "Node.h"

using namespace std;

//! �Զ���ElementTypes���������ͣ���ö�ٳ����п��ܳ��ֵĵ�Ԫ����
enum ElementTypes
{
    UNDEFINED = 0,
    Bar,    // Bar element
    Q4,     // 4Q element
    T3,     // 3T element
    H8,     // 8H element
    Beam,   // Beam element
    Plate,  // Plate element
    Shell   // Shell elment
};

//! ��Ԫ����(һ����Ԫ������ֻ��һ�ֵ�Ԫ)
class CElementGroup
{
private:
    //! �������нڵ���б�(�������������Ԫ���漰�Ľڵ�)������������κ�һ����Ԫ�飬��һ����̬����
    static CNode* NodeList_;

    //! ����õ�Ԫ��ĵ�Ԫ����
    ElementTypes ElementType_;

    //! �õ�Ԫ���һ����Ԫ��ռ���ڴ��С
    std::size_t ElementSize_;

    //! �õ�Ԫ��ĵ�Ԫ����
    unsigned int NUME_;

    //! ����һ��CElement���͵�ָ��
    CElement* ElementList_;

    //! �����������
    unsigned int NUMMAT_;

    //! ����һ��CMaterial���͵�ָ��
    CMaterial* MaterialList_;

    //! �õ�Ԫ���ĳ�ֲ�����ռ���ڴ��С(����ϱ����޹أ�ֻ�뵥Ԫ�����й�)
    std::size_t MaterialSize_;

public:
    //! ���캯��
    CElementGroup();

    //! ��������
    ~CElementGroup();

    //! 
    bool Read(ifstream& Input);

    //! Elementtype_�����󣬸ú����������ElementSize_��MaterialSize_��Ϊ���������[]��GetMaterial������ʵ�ַ���
    void CalculateMemberSize();

    //! ���뵥Ԫ����������һ������������
    void AllocateElements(std::size_t size);

    //! �����������������һ������������
    void AllocateMaterials(std::size_t size);

    //! ���ص�Ԫ���Ӧ�ĵ�Ԫ����
    ElementTypes GetElementType() { return ElementType_; }

    //! ���ص�Ԫ���еĵ�Ԫ����
    unsigned int GetNUME() { return NUME_; }

    //������������������ڴ�ElementGroup�����ȡ��i��Element����
    CElement& operator[](unsigned int i);

    //! ���ڴ�ElementGroup�����ȡ��i�����϶���
    CMaterial& GetMaterial(unsigned int i);

    //! ���ص�Ԫ���в�������
    unsigned int GetNUMMAT() { return NUMMAT_; }
};
