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

//��̬������ʼ��
CNode* CElementGroup::NodeList_ = nullptr;

//! ���캯��
CElementGroup::CElementGroup()
{
	//��CDomain�����л�ȡNodeList
    if (!NodeList_)
    {
        CDomain* FEMData = CDomain::GetInstance();
        NodeList_ = FEMData->GetNodeList();
    }
    //��ʼ����Ԫ���ͣ���Ԫ��������Ԫ�б����������������б�
    ElementType_ = ElementTypes::UNDEFINED;
    
    NUME_ = 0;
    ElementList_ = nullptr;
    
    NUMMAT_ = 0;
    MaterialList_ = nullptr;
}

//! ���������Ȳ��ù�
CElementGroup::~CElementGroup()
{
    if (ElementList_)
        delete [] ElementList_;
    
    if (MaterialList_)
        delete [] MaterialList_;
}

//! ��[]��������������ء�
//  ���뵥Ԫ���еĵ�Ԫ���i,����õ�Ԫ��������ʹ�ø�������ý�������Ķ��󣬻ᷢ����̬

CElement& CElementGroup::operator[](unsigned int i)
{
    return *(CElement*)((std::size_t)(ElementList_) + i*ElementSize_);
	//  �﷨������*(CElement*)((std::size_t)(ElementList_) + i*ElementSize_);
//��һ��*  �����ã���ָ��CElement�����ָ����ָ��ָ��Ķ���
//(CElement*)  ǿ������ת��������������ת����ָ������
//(std::size_t)   ǿ������ת������ָ������ת�����������ͣ������������
}

//  ����GetMaterial���������뵥Ԫ���еĲ��ϱ��i��������϶��󡣻ᷢ����̬
CMaterial& CElementGroup::GetMaterial(unsigned int i)
{
    return *(CMaterial*)((std::size_t)(MaterialList_) + i*MaterialSize_);
}

//! ���㵥Ԫ��ĵ�Ԫ�͸õ�Ԫ��Ӧ�Ĳ�����ռ���ڴ��С�����ѽ���洢��CElementGroup������ElementSize_��MaterialSize_�С�
//  ע�������ռ�ڴ�Ĵ�С��ʲô�����޹أ���ȫȡ���ڵ�Ԫ����
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

//! ���뵥Ԫ�����ڶ�������һ��CBar(�����������͵�Ԫ)���͵����飬�洢��ElementList_��
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

//! ��������������ڶ�������һ��CBarmaterial(�����������͵�Ԫ)���͵����飬�浽MaterialList_��
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

//! ����һ��ifstream�����ṩһ���ļ��ӿڣ����ļ��ж�ȡ����
bool CElementGroup::Read(ifstream& Input)
{
	//��ȡ��Ԫ���͡���Ԫ���е�Ԫ������Ԫ���в�����
    Input >> (int&)ElementType_ >> NUME_ >> NUMMAT_;
    //����ElementGroup���е�CalculateMemberSize()���������ݵ�Ԫ���ͼ��㵥Ԫ��������ռ���ڴ�
    CalculateMemberSize();

//  �ڶ�������һ������ΪCMaterial����ΪNUMMAT_������
    AllocateMaterials(NUMMAT_);
    
//  ��ȡ���еĲ��ϲ���
    for (unsigned int mset = 0; mset < NUMMAT_; mset++)
    {
		//GetMaterial������ĸ�������ã�����return����������󣬷�����̬�������Read��д�������е�Read��
        GetMaterial(mset).Read(Input);
        //�ж������ʽ�Ƿ���ȷ
        if (GetMaterial(mset).nset != mset + 1)
        {
            cerr << "*** Error *** Material sets must be inputted in order !" << endl
            << "    Expected set : " << mset + 1 << endl
            << "    Provided set : " << GetMaterial(mset).nset << endl;
        
            return false;
        }
    }

//  �ڶ�������һ������ΪCElement����ΪNUME_������
    AllocateElements(NUME_);
    
//  ���뵥Ԫ����
    for (unsigned int Ele = 0; Ele < NUME_; Ele++)
    {
        unsigned int N;
        
        Input >> N;    // element number
        //�жϵ�Ԫ����Ƿ���˳����룬������Ǿͱ���
        if (N != Ele + 1)
        {
            cerr << "*** Error *** Elements must be inputted in order !" << endl
            << "    Expected element : " << Ele + 1 << endl
            << "    Provided element : " << N << endl;
            
            return false;
        }
		//ʹ�������ص�[]�������
        if (!(*this)[Ele].Read(Input, MaterialList_, NodeList_))
            return false;
    }

    return true;
}
