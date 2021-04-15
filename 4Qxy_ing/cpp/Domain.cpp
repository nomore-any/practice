/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Domain.h"
#include "Material.h"

using namespace std;

//	清除一个数组
template <class type> void clear( type* a, unsigned int N )
{
	for (unsigned int i = 0; i < N; i++)
		a[i] = 0;
}

CDomain* CDomain::_instance = nullptr;

//	构造函数，初始化
CDomain::CDomain()
{
	Title[0] = '0';
	MODEX = 0;

	NUMNP = 0;
	NodeList = nullptr;
	
	NUMEG = 0;
	EleGrpList = nullptr;
	
	NLCASE = 0;
	NLOAD = nullptr;
	LoadCases = nullptr;
	
	NEQ = 0;

	Force = nullptr;
	StiffnessMatrix = nullptr;
}

//	析构函数
CDomain::~CDomain()
{
	delete [] NodeList;

	delete [] EleGrpList;

	delete [] NLOAD;
	delete [] LoadCases;

	delete [] Force;
	delete StiffnessMatrix;
}

//	如果_instance是一个空指针，说明还没有CDomain对象，此时在堆区生成一个CDomain对象，返回指向该对象的指针
CDomain* CDomain::GetInstance()
{
	if (!_instance) 
		_instance = new CDomain();
	
	return _instance;
}

//	第一个传入的参数是希望读取的文件名，第二个传入的参数是希望被写入的文件名(outputter实现两个目标，一个是直接打印，一个是把结果写到一个文件里)
bool CDomain::ReadData(string FileName, string OutFile)
{
	//用Input读取文件
	Input.open(FileName);
	//判断文件是否存在，不用管
	if (!Input) 
	{
		cerr << "*** Error *** File " << FileName << " does not exist !" << endl;
		exit(3);
	}
	//使用COutputter类里的GetInstance函数，生成一个Outputter对象用于输出所有的东西，并产生写入文件的接口，返回指针
	COutputter* Output = COutputter::GetInstance(OutFile);

//	读取标题
	Input.getline(Title, 256);
	Output->OutputHeading();

//	读取如下各种全局信息
	Input >> NUMNP >> NUMEG >> NLCASE >> MODEX;

//	用CDomain类的ReadNodalPoints函数读取节点信息，并且打印出来
	if (ReadNodalPoints())
        Output->OutputNodeInfo();
    else
        return false;

//	用CDomain类的函数计算NEQ，用Outputter类输出NEQ
	CalculateEquationNumber();
	Output->OutputEquationNumber();

//	用CDomain类的函数读取载荷信息，用Outputter类输出
	if (ReadLoadCases())
        Output->OutputLoadInfo();
    else
        return false;

//	用CDomain类的函数读取单元信息，用Outputter类输出
	if (ReadElements())
        Output->OutputElementInfo();
    else
        return false;

	return true;
}

//	读取节点信息
bool CDomain::ReadNodalPoints()
{

//	在堆区产生一个长度为'节点总数'的CNode类型数组
	NodeList = new CNode[NUMNP];

//	对所有的节点循环
	for (unsigned int np = 0; np < NUMNP; np++)
    {
		//这一步用CNode类中的Read函数读取数据，既起到读取的作用，也起到判断输入格式是否正确的作用
		if (!NodeList[np].Read(Input))
			return false;
        //判断节点是否按照节点编号的顺序输入
        if (NodeList[np].NodeNumber != np + 1)
        {
            cerr << "*** Error *** Nodes must be inputted in order !" << endl
            << "   Expected node number : " << np + 1 << endl
            << "   Provided node number : " << NodeList[np].NodeNumber << endl;
        
            return false;
        }
    }

	return true;
}

//	这个函数起到两个作用，一是计算总方程数，二是把bcode的值变成自由度全局编号
void CDomain::CalculateEquationNumber()
{
	NEQ = 0;
	for (unsigned int np = 0; np < NUMNP; np++)	// Loop over for all node
	{
		for (unsigned int dof = 0; dof < CNode::NDF; dof++)	// Loop over for DOFs of node np
		{
			if (NodeList[np].bcode[dof]) 
				NodeList[np].bcode[dof] = 0;
			else
			{
				NEQ++;
				NodeList[np].bcode[dof] = NEQ;
			}
		}
	}
}

//	读取载荷信息
bool CDomain::ReadLoadCases()
{
//	分配载荷工况数的CLoadCaseData类型数组
	LoadCases = new CLoadCaseData[NLCASE];	// List all load cases

//	对载荷工况编号循环，调用CLoadCaseData类中的Read函数读取该工况的载荷信息
	for (unsigned int lcase = 0; lcase < NLCASE; lcase++)
    {
        unsigned int LL;
        Input >> LL;
        
		//如果不按编号顺序输入载荷会报错
        if (LL != lcase + 1)
        {
            cerr << "*** Error *** Load case must be inputted in order !" << endl
            << "   Expected load case : " << lcase + 1 << endl
            << "   Provided load case : " << LL << endl;
            
            return false;
        }

        LoadCases[lcase].Read(Input);
    }

	return true;
}

// 读取单元信息
bool CDomain::ReadElements()
{
	//分配单元组类型数组
    EleGrpList = new CElementGroup[NUMEG];

//	对单元组编号循环
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
		//用ElementGroup类中的Read函数读取单元数据，同时判断输入格式是否正确
        if (!EleGrpList[EleGrp].Read(Input))
            return false;
    
    return true;
}

//	计算全局刚度阵的列高
void CDomain::CalculateColumnHeights()
{
#ifdef _DEBUG_
    COutputter* Output = COutputter::GetInstance();
    *Output << setw(9) << "Ele = " << setw(22) << "Location Matrix" << endl;
#endif

	//思路是，从ElementGroup开始循环，再对Element循环，把Element的参数传进CSkylineMatrix的CalculateColumnHeight函数，求得列高数组
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)		//	Loop over for all element groups
    {
        CElementGroup& ElementGrp = EleGrpList[EleGrp];
        unsigned int NUME = ElementGrp.GetNUME();
        
		for (unsigned int Ele = 0; Ele < NUME; Ele++)	//	Loop over for all elements in group EleGrp
        {
            CElement& Element = ElementGrp[Ele];

            // 在这一步给每个单元生成分配数组(LocationMatrix，局部自由度和全局自由度的编号联系)
            Element.GenerateLocationMatrix();
            
#ifdef _DEBUG_
            unsigned int* LocationMatrix = Element.GetLocationMatrix();
            
            *Output << setw(9) << Ele+1;
            for (int i=0; i<Element.GetND(); i++)
                *Output << setw(5) << LocationMatrix[i];
            *Output << endl;
#endif
			//计算列高，存进StiffnessMatrix的ColumnHeights_变量里
            StiffnessMatrix->CalculateColumnHeight(Element.GetLocationMatrix(), Element.GetND());
        }
    }
    //计算最大半带宽
    StiffnessMatrix->CalculateMaximumHalfBandwidth();
    
#ifdef _DEBUG_
    *Output << endl;
	Output->PrintColumnHeights();
#endif

}

//    为各种矩阵和向量分配内存空间，并且调用CalculateColumnHeight函数获取列高，输出一些与刚度阵有关的信息
void CDomain::AllocateMatrices()
{
    //    Allocate for global force/displacement vector
    Force = new double[NEQ];
    
    //  Create the banded stiffness matrix
    StiffnessMatrix = new CSkylineMatrix<double>(NEQ);
    
    //    Calculate column heights
    CalculateColumnHeights();
    
    //    Calculate address of diagonal elements in banded matrix
    StiffnessMatrix->CalculateDiagnoalAddress();
    
    //    Allocate for banded global stiffness matrix
    StiffnessMatrix->Allocate();
    
    COutputter* Output = COutputter::GetInstance();
    Output->OutputTotalSystemData();
}

//	组装全局刚度阵
void CDomain::AssembleStiffnessMatrix()
{
//	思路是，对单元组循环，在每个单元组内对单元循环，对每个单元调用ElementStiffness函数算得单元刚度阵，传入StiffnessMatrix的Assembly函数里，实现组装
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
		//单元组循环
        CElementGroup& ElementGrp = EleGrpList[EleGrp];
        unsigned int NUME = ElementGrp.GetNUME();

		unsigned int size = ElementGrp[0].SizeOfStiffnessMatrix();
		double* Matrix = new double[size];

//		单元循环
		for (unsigned int Ele = 0; Ele < NUME; Ele++)
        {
            CElement& Element = ElementGrp[Ele];
			//计算单元刚度阵存到Matrix矩阵，组装进全局刚度阵
            Element.ElementStiffness(Matrix);
            StiffnessMatrix->Assembly(Matrix, Element.GetLocationMatrix(), Element.GetND());
        }

		//清除堆区内存，先不管
		delete[] Matrix;
		Matrix = nullptr;
	}

#ifdef _DEBUG_
	COutputter* Output = COutputter::GetInstance();
	Output->PrintStiffnessMatrix();
#endif

}

//	传入一个工况编号，生成该工况的全局力向量
bool CDomain::AssembleForce(unsigned int LoadCase)
{
	if (LoadCase > NLCASE) 
		return false;

	CLoadCaseData* LoadData = &LoadCases[LoadCase - 1];

	//不同工况是不叠加的，每个载荷工况都会单独求解位移，所以这里把Force清零，再计算下一个工况
    clear(Force, NEQ);

//	对该载荷工况的载荷编号进行循环
	for (unsigned int lnum = 0; lnum < LoadData->nloads; lnum++)
	{
		//找到第lnum个载荷所加自由度的全局编号
		//方法是Domain提供NodeList,LoadData提供节点编号，找到节点，然后用节点的bcode找到自由度的全局编号
		unsigned int dof = NodeList[LoadData->node[lnum] - 1].bcode[LoadData->dof[lnum] - 1];
        
        if(dof) // The DOF is activated
            Force[dof - 1] += LoadData->load[lnum];
	}

	return true;
}

