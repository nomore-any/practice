/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "ElementGroup.h"
#include "Domain.h"

//静态变量初始化
CNode* CElementGroup::NodeList_ = nullptr;

//! 构造函数
CElementGroup::CElementGroup()
{
	//从CDomain对象中获取NodeList
    if (!NodeList_)
    {
        CDomain* FEMData = CDomain::GetInstance();
        NodeList_ = FEMData->GetNodeList();
    }
    //初始化单元类型，单元总数，单元列表，材料总数，材料列表
    ElementType_ = ElementTypes::UNDEFINED;
    
    NUME_ = 0;
    ElementList_ = nullptr;
    
    NUMMAT_ = 0;
    MaterialList_ = nullptr;
}

//! 析构函数先不用管
CElementGroup::~CElementGroup()
{
    if (ElementList_)
        delete [] ElementList_;
    
    if (MaterialList_)
        delete [] MaterialList_;
}

//! 对[]进行了运算符重载。
//  输入单元组中的单元编号i,输出该单元对象，由于使用父类的引用接收子类的对象，会发生多态

CElement& CElementGroup::operator[](unsigned int i)
{
    return *(CElement*)((std::size_t)(ElementList_) + i*ElementSize_);
	//  语法分析：*(CElement*)((std::size_t)(ElementList_) + i*ElementSize_);
//第一个*  解引用，把指向CElement对象的指针变成指针指向的对象
//(CElement*)  强制类型转换，把整数类型转换成指针类型
//(std::size_t)   强制类型转换，把指针类型转换成整数类型，便于求和运算
}

//  定义GetMaterial函数，输入单元组中的材料编号i，输出材料对象。会发生多态
CMaterial& CElementGroup::GetMaterial(unsigned int i)
{
    return *(CMaterial*)((std::size_t)(MaterialList_) + i*MaterialSize_);
}

//! 计算单元组的单元和该单元对应的材料所占的内存大小，并把结果存储到CElementGroup的属性ElementSize_、MaterialSize_中。
//  注意材料所占内存的大小与什么材料无关，完全取决于单元类型
void CElementGroup::CalculateMemberSize()
{
    switch (ElementType_)
    {
        case ElementTypes::UNDEFINED:
            std::cerr << "Setting element type to UNDEFINED." << std::endl;
            exit(5);
        case ElementTypes::Bar:
            ElementSize_ = sizeof(CBar);
            MaterialSize_ = sizeof(CBarMaterial);
            break;
        case ElementTypes::Q4:
            ElementSize_ = sizeof(CQ4);
            MaterialSize_ = sizeof(CQ4Material);
            break;
        default:
            std::cerr << "Type " << ElementType_ << " not available. See CElementGroup::CalculateMemberSize." << std::endl;
            exit(5);
            break;
    }
}

//! 传入单元数，在堆区分配一个CBar(或者其他类型单元)类型的数组，存储到ElementList_中
void CElementGroup::AllocateElements(std::size_t size)
{
    switch(ElementType_)
    {
        case ElementTypes::Bar:
            ElementList_ = new CBar[size];
            break;
        default:
            std::cerr << "Type " << ElementType_ << " not available. See CElementGroup::AllocateElement." << std::endl;
            exit(5);
    }
}

//! 传入材料总数，在堆区分配一个CBarmaterial(或者其他类型单元)类型的数组，存到MaterialList_中
void CElementGroup::AllocateMaterials(std::size_t size)
{
    switch(ElementType_)
    {
        case ElementTypes::Bar:
            MaterialList_ = new CBarMaterial[size];
            break;
        default:
            std::cerr << "Type " << ElementType_ << " not available. See CElementGroup::AllocateMaterial." << std::endl;
            exit(5);
    }
}

//! 传入一个ifstream对象，提供一个文件接口，从文件中读取数据
bool CElementGroup::Read(ifstream& Input)
{
	//读取单元类型、单元组中单元数、单元组中材料数
    Input >> (int&)ElementType_ >> NUME_ >> NUMMAT_;
    //调用ElementGroup类中的CalculateMemberSize()函数，根据单元类型计算单元、材料所占的内存
    CalculateMemberSize();

//  在堆区产生一个类型为CMaterial长度为NUMMAT_的数组
    AllocateMaterials(NUMMAT_);
    
//  读取所有的材料参数
    for (unsigned int mset = 0; mset < NUMMAT_; mset++)
    {
		//GetMaterial的输出的父类的引用，但是return的是子类对象，发生多态。下面的Read是写在子类中的Read。
        GetMaterial(mset).Read(Input);
        //判断输入格式是否正确
        if (GetMaterial(mset).nset != mset + 1)
        {
            cerr << "*** Error *** Material sets must be inputted in order !" << endl
            << "    Expected set : " << mset + 1 << endl
            << "    Provided set : " << GetMaterial(mset).nset << endl;
        
            return false;
        }
    }

//  在堆区产生一个类型为CElement长度为NUME_的数组
    AllocateElements(NUME_);
    
//  读入单元数据
    for (unsigned int Ele = 0; Ele < NUME_; Ele++)
    {
        unsigned int N;
        
        Input >> N;    // element number
        //判断单元编号是否按照顺序读入，如果不是就报错
        if (N != Ele + 1)
        {
            cerr << "*** Error *** Elements must be inputted in order !" << endl
            << "    Expected element : " << Ele + 1 << endl
            << "    Provided element : " << N << endl;
            
            return false;
        }
		//使用了重载的[]运算符，
        if (!(*this)[Ele].Read(Input, MaterialList_, NodeList_))
            return false;
    }

    return true;
}
