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

//! 自定义ElementTypes的数据类型，并枚举出所有可能出现的单元种类
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

//! 单元组类(一个单元组里面只有一种单元)
class CElementGroup
{
private:
    //! 包含所有节点的列表(不仅仅是这个单元组涉及的节点)，因而不属于任何一个单元组，是一个静态变量
    static CNode* NodeList_;

    //! 定义该单元组的单元类型
    ElementTypes ElementType_;

    //! 该单元组的一个单元所占的内存大小
    std::size_t ElementSize_;

    //! 该单元组的单元总数
    unsigned int NUME_;

    //! 定义一个CElement类型的指针
    CElement* ElementList_;

    //! 定义材料总数
    unsigned int NUMMAT_;

    //! 定义一个CMaterial类型的指针
    CMaterial* MaterialList_;

    //! 该单元组的某种材料所占的内存大小(与材料本身无关，只与单元类型有关)
    std::size_t MaterialSize_;

public:
    //! 构造函数
    CElementGroup();

    //! 析构函数
    ~CElementGroup();

    //! 
    bool Read(ifstream& Input);

    //! Elementtype_给定后，该函数用于求出ElementSize_与MaterialSize_，为重载运算符[]、GetMaterial函数的实现服务
    void CalculateMemberSize();

    //! 传入单元总数，申请一个堆区的数组
    void AllocateElements(std::size_t size);

    //! 传入材料总数，申请一个堆区的数组
    void AllocateMaterials(std::size_t size);

    //! 返回单元组对应的单元类型
    ElementTypes GetElementType() { return ElementType_; }

    //! 返回单元组中的单元个数
    unsigned int GetNUME() { return NUME_; }

    //声明重载运算符，用于从ElementGroup对象获取第i个Element对象
    CElement& operator[](unsigned int i);

    //! 用于从ElementGroup对象获取第i个材料对象
    CMaterial& GetMaterial(unsigned int i);

    //! 返回单元组中材料总数
    unsigned int GetNUMMAT() { return NUMMAT_; }
};
