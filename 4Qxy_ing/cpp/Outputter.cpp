/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include <ctime>

#include "Domain.h"
#include "Outputter.h"
#include "SkylineMatrix.h"

using namespace std;

//	打印时间，不用管
void COutputter::PrintTime(const struct tm* ptm, COutputter &output)
{
	const char* weekday[] = {"Sunday", "Monday", "Tuesday", "Wednesday",
							 "Thursday", "Friday", "Saturday"};
	const char* month[] = {"January", "February", "March", "April", "May", "June",
						   "July", "August", "September", "October", "November", "December"};

	output << "        (";
	output << ptm->tm_hour << ":" << ptm->tm_min << ":" << ptm->tm_sec << " on ";
	output << month[ptm->tm_mon] << " " << ptm->tm_mday << ", " << ptm->tm_year + 1900 << ", "
		   << weekday[ptm->tm_wday] << ")" << endl
		   << endl;
}
//对COutputter类的属性_instance进行初始化，因为这个属性不是属于任何一个对象的，所以不能写在构造函数里，而是单独初始化
COutputter* COutputter::_instance = nullptr;

//	输入文件名，用ofstream对象打开文件
COutputter::COutputter(string FileName)
{
	OutputFile.open(FileName);

	if (!OutputFile)
	{
		cerr << "*** Error *** File " << FileName << " does not exist !" << endl;
		exit(3);
	}
}

//	当_instance是空指针的时候，让它指向一个新建的堆区的COutputter对象，并返回该指针
//  如果希望其他类的函数调用Outputter对象，有两种办法，一种就是使用该指针，另一种就是直接传入outputter对象，基本上都是后者
COutputter* COutputter::GetInstance(string FileName)
{
	if (!_instance)
		_instance = new COutputter(FileName);
    
	return _instance;
}

//	不用管
void COutputter::OutputHeading()
{
	CDomain* FEMData = CDomain::GetInstance();

	*this << "TITLE : " << FEMData->GetTitle() << endl;

	time_t rawtime;
	struct tm* timeinfo;

	time(&rawtime);
	timeinfo = localtime(&rawtime);

	PrintTime(timeinfo, *this);
}

//	打印节点信息
void COutputter::OutputNodeInfo()
{
	//从CDomain获取指向CDomain对象的指针
	CDomain* FEMData = CDomain::GetInstance();
	//从CDomain获取CNode类型的数组NodeList，即节点数组
	CNode* NodeList = FEMData->GetNodeList();

	//不用管
	*this << "C O N T R O L   I N F O R M A T I O N" << endl
		  << endl;
	*this << setiosflags(ios::scientific) << setprecision(5);

	//利用FEMData指针获取CDomain对象中的'节点总数'、'单元组个数'、'载荷工况数'、'求解模式'
	unsigned int NUMNP = FEMData->GetNUMNP();
	unsigned int NUMEG = FEMData->GetNUMEG();
	unsigned int NLCASE = FEMData->GetNLCASE();
	unsigned int MODEX = FEMData->GetMODEX();

	//打印'节点总数'、'单元组个数'、'载荷工况数'、'求解模式'
	*this << "      NUMBER OF NODAL POINTS . . . . . . . . . . (NUMNP)  =" << setw(6) << NUMNP << endl;
	*this << "      NUMBER OF ELEMENT GROUPS . . . . . . . . . (NUMEG)  =" << setw(6) << NUMEG << endl;
	*this << "      NUMBER OF LOAD CASES . . . . . . . . . . . (NLCASE) =" << setw(6) << NLCASE << endl;
	*this << "      SOLUTION MODE  . . . . . . . . . . . . . . (MODEX)  =" << setw(6) << MODEX << endl;
	*this << "         EQ.0, DATA CHECK" << endl
		  << "         EQ.1, EXECUTION" << endl
		  << endl;

	*this << " N O D A L   P O I N T   D A T A" << endl << endl;
	*this << "    NODE       BOUNDARY                         NODAL POINT" << endl
		  << "   NUMBER  CONDITION  CODES                     COORDINATES" << endl;

	//利用CNode类中的Write函数把各种数据打印出来以及写进目标文件里
	for (unsigned int np = 0; np < NUMNP; np++)
		NodeList[np].Write(*this);

	*this << endl;
}

//	对节点编号进行循环，输出每个节点的3个自由度对应的全局编号(EquationNumber)
void COutputter::OutputEquationNumber()
{
	//从CDomain获取指向CDomain对象的指针
	CDomain* FEMData = CDomain::GetInstance();
	//使用FEMData指针获取'节点总数'
	unsigned int NUMNP = FEMData->GetNUMNP();
	//使用FEMData指针获取CDomain对象中CNode类型的数组NodeList
	CNode* NodeList = FEMData->GetNodeList();

	*this << " EQUATION NUMBERS" << endl
		  << endl;
	*this << "   NODE NUMBER   DEGREES OF FREEDOM" << endl;
	*this << "        N           X    Y    Z" << endl;
	//对节点循环，利用CNode类中的WriteEquationNo输出该节点三个自由度的全局编号
	for (unsigned int np = 0; np < NUMNP; np++) // Loop over for all node
		NodeList[np].WriteEquationNo(*this);

	*this << endl;
}

//	输出单元信息
void COutputter::OutputElementInfo()
{
	//获取指向CDomain对象的指针
	CDomain* FEMData = CDomain::GetInstance();
	//从CDomain对象中获取单元组数
	unsigned int NUMEG = FEMData->GetNUMEG();

	*this << " E L E M E N T   G R O U P   D A T A" << endl
		  << endl
		  << endl;
	//对单元组编号进行循环
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
		*this << " E L E M E N T   D E F I N I T I O N" << endl
			  << endl;
		//为了获取该单元组的单元类型，先从单元组数组得到单元组，再用单元组类中的GetElementType()函数获取单元类型
		//同理也用这种方法获取单元组中的单元总数NUME
		ElementTypes ElementType = FEMData->GetEleGrpList()[EleGrp].GetElementType();
		unsigned int NUME = FEMData->GetEleGrpList()[EleGrp].GetNUME();

		*this << " ELEMENT TYPE  . . . . . . . . . . . . .( NPAR(1) ) . . =" << setw(5)
			  << ElementType << endl;
		*this << "     EQ.1, TRUSS ELEMENTS" << endl
			  << "     EQ.2, ELEMENTS CURRENTLY" << endl
			  << "     EQ.3, NOT AVAILABLE" << endl
			  << endl;

		*this << " NUMBER OF ELEMENTS. . . . . . . . . . .( NPAR(2) ) . . =" << setw(5) << NUME
			  << endl
			  << endl;

		switch (ElementType)
		{
			case ElementTypes::Bar: // Bar element
				//函数定义在下面。输出该ElementGroup中所有单元的信息
				OutputBarElements(EleGrp);
				break;
			case ElementTypes::Q4: // Q4 element
				//函数定义在下面。输出该ElementGroup中所有单元的信息
				OutputQ4Elements(EleGrp);
				break;
		    default:
		        *this << ElementType << " has not been implemented yet." << endl;
		        break;
		}
	}
}
//    传入单元组编号(必须是杆单元的单元组)，输出该单元组所有单元的信息
void COutputter::OutputBarElements(unsigned int EleGrp)
{
	//获取指向CDomain对象的指针
	CDomain* FEMData = CDomain::GetInstance();
	//获取该单元组对象以及单元组的材料总数
	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();
	//输出材料信息
	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		  << endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
	*this << " AND CROSS-SECTIONAL  CONSTANTS  . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
		  << endl
		  << endl;

	*this << "  SET       YOUNG'S     CROSS-SECTIONAL" << endl
		  << " NUMBER     MODULUS          AREA" << endl
		  << "               E              A" << endl;
	//没看懂@
	*this << setiosflags(ios::scientific) << setprecision(5);

	for (unsigned int mset = 0; mset < NUMMAT; mset++)
    {
        *this << setw(5) << mset+1;
	//输出材料参数的方法：使用ElementGroup类中的GetMaterial函数找到单元组中某种材料的对象，利用CBarMaterial类中的Write函数输出参数
		ElementGroup.GetMaterial(mset).Write(*this);
    }

	*this << endl << endl
		  << " E L E M E N T   I N F O R M A T I O N" << endl;
    
	*this << " ELEMENT     NODE     NODE       MATERIAL" << endl
		  << " NUMBER-N      I        J       SET NUMBER" << endl;
	//输出所有单元信息，原理同上
	unsigned int NUME = ElementGroup.GetNUME();

	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < NUME; Ele++)
    {
        *this << setw(5) << Ele+1;
		ElementGroup[Ele].Write(*this);
    }

	*this << endl;
}

//    传入单元组编号(必须是4Q单元的单元组)，输出该单元组所有单元的信息
void COutputter::OutputQ4Elements(unsigned int EleGrp)
{
	//获取指向CDomain对象的指针
	CDomain* FEMData = CDomain::GetInstance();
	//获取该单元组对象以及单元组的材料总数
	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();
	//输出材料信息
	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		  << endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
	*this << " AND CROSS-SECTIONAL  CONSTANTS  . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
		  << endl
		  << endl;

	*this << "  SET       YOUNG'S     CROSS-SECTIONAL" << endl
		  << " NUMBER     MODULUS          泊松比" << endl
		  << "               E              Posion" << endl;
	//没看懂@
	*this << setiosflags(ios::scientific) << setprecision(5);

	for (unsigned int mset = 0; mset < NUMMAT; mset++)
    {
        *this << setw(5) << mset+1;
	//输出材料参数的方法：使用ElementGroup类中的GetMaterial函数找到单元组中某种材料的对象，利用CQ4Material类中的Write函数输出参数
		ElementGroup.GetMaterial(mset).Write(*this);
    }

	*this << endl << endl
		  << " E L E M E N T   I N F O R M A T I O N" << endl;
    
	*this << " ELEMENT     NODE     NODE     NODE     NODE     MATERIAL" << endl
		  << " NUMBER-N     I        J        K        L      SET NUMBER" << endl;
	//输出所有单元信息，原理同上
	unsigned int NUME = ElementGroup.GetNUME();

	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < NUME; Ele++)
    {
        *this << setw(5) << Ele+1;
		ElementGroup[Ele].Write(*this);
    }

	*this << endl;
}

//	打印载荷信息
void COutputter::OutputLoadInfo()
{
	//获取指向CDomain对象的指针
	CDomain* FEMData = CDomain::GetInstance();
	//对工况编号循环
	for (unsigned int lcase = 1; lcase <= FEMData->GetNLCASE(); lcase++)
	{
		//获取指向某个CLoadCaseData对象的指针
		CLoadCaseData* LoadData = &FEMData->GetLoadCases()[lcase - 1];

		//输出各种东西
		*this << setiosflags(ios::scientific);
		*this << " L O A D   C A S E   D A T A" << endl
			  << endl;

		*this << "     LOAD CASE NUMBER . . . . . . . =" << setw(6) << lcase << endl;
		*this << "     NUMBER OF CONCENTRATED LOADS . =" << setw(6) << LoadData->nloads << endl
			  << endl;
		*this << "    NODE       DIRECTION      LOAD" << endl
			  << "   NUMBER                   MAGNITUDE" << endl;
		//使用CLoadCaseData类下的Write函数输出
		LoadData->Write(*this);

		*this << endl;
	}
}

//	输出位移信息
void COutputter::OutputNodalDisplacement()
{
	//从CDomain类获取指针、节点列表、位移数组
	CDomain* FEMData = CDomain::GetInstance();
	CNode* NodeList = FEMData->GetNodeList();
	double* Displacement = FEMData->GetDisplacement();

	*this << setiosflags(ios::scientific);

	*this << " D I S P L A C E M E N T S" << endl
		  << endl;
	*this << "  NODE           X-DISPLACEMENT    Y-DISPLACEMENT    Z-DISPLACEMENT" << endl;
	//输出每个节点的位移信息，使用CNode类下的WriteNodalDisplacement函数
	for (unsigned int np = 0; np < FEMData->GetNUMNP(); np++)
		NodeList[np].WriteNodalDisplacement(*this, Displacement);

	*this << endl;
}

//	输出单元应力
void COutputter::OutputElementStress()
{
	//获取指针
	CDomain* FEMData = CDomain::GetInstance();
	//获取全局位移向量
	double* Displacement = FEMData->GetDisplacement();
	//获取单元组个数
	unsigned int NUMEG = FEMData->GetNUMEG();
    //对单元组循环
	for (unsigned int EleGrpIndex = 0; EleGrpIndex < NUMEG; EleGrpIndex++)
	{
		*this << " S T R E S S  C A L C U L A T I O N S  F O R  E L E M E N T  G R O U P" << setw(5)
			  << EleGrpIndex + 1 << endl
			  << endl;
		//从CDomain类获取单元组
		CElementGroup& EleGrp = FEMData->GetEleGrpList()[EleGrpIndex];
		//从单元组获取单元总数、单元类型
		unsigned int NUME = EleGrp.GetNUME();
		ElementTypes ElementType = EleGrp.GetElementType();

		switch (ElementType)
		{
			case ElementTypes::Bar: // Bar element
				*this << "  ELEMENT             FORCE            STRESS" << endl
					<< "  NUMBER" << endl;

				double stress;
				//对单元循环，输出每个单元的应力
				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					CElement& Element = EleGrp[Ele];
					Element.ElementStress(&stress, Displacement);
		//由于GetElementMaterial()是写在父类里面的，所以指向CMaterial类，需要进行强制类型转换，把结果变为指向CBarMaterial类
					CBarMaterial& material = *dynamic_cast<CBarMaterial*>(Element.GetElementMaterial());
					*this << setw(5) << Ele + 1 << setw(22) << stress * material.Area << setw(18)
						<< stress << endl;
				}

				*this << endl;

				break;

			case ElementTypes::Q4: // 4Q element
				*this << "  ELEMENT             FORCE            STRESS" << endl
					<< "  NUMBER" << endl;

				double stress;
				//对单元循环，输出每个单元的应力
				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					CElement& Element = EleGrp[Ele];
					Element.ElementStress(&stress, Displacement);
		//由于GetElementMaterial()是写在父类里面的，所以指向CMaterial类，需要进行强制类型转换，把结果变为指向CBarMaterial类
					CQ4Material& material = *dynamic_cast<CQ4Material*>(Element.GetElementMaterial());
					*this << setw(5) << Ele + 1 << setw(22) << stress * material. << setw(18)
						<< stress << endl;
				}

				*this << endl;

				break;

			default: // Invalid element type
				cerr << "*** Error *** Elment type " << ElementType
					<< " has not been implemented.\n\n";
		}
	}
}

//	输出一些与刚度阵有关的东西
void COutputter::OutputTotalSystemData()
{
	CDomain* FEMData = CDomain::GetInstance();

	*this << "	TOTAL SYSTEM DATA" << endl
		  << endl;

	*this << "     NUMBER OF EQUATIONS . . . . . . . . . . . . . .(NEQ) = " << FEMData->GetNEQ()
		  << endl
		  << "     NUMBER OF MATRIX ELEMENTS . . . . . . . . . . .(NWK) = " << FEMData->GetStiffnessMatrix()->size()
		  << endl
		  << "     MAXIMUM HALF BANDWIDTH  . . . . . . . . . . . .(MK ) = " << FEMData->GetStiffnessMatrix()->GetMaximumHalfBandwidth()
		  << endl
		  << "     MEAN HALF BANDWIDTH . . . . . . . . . . . . . .(MM ) = " << FEMData->GetStiffnessMatrix()->size() / FEMData->GetNEQ() << endl
		  << endl
		  << endl;
}

#ifdef _DEBUG_

//	Print column heights for debuging
void COutputter::PrintColumnHeights()
{
	*this << "*** _Debug_ *** Column Heights" << endl;

	CDomain* FEMData = CDomain::GetInstance();

	unsigned int NEQ = FEMData->GetNEQ();
	CSkylineMatrix<double> *StiffnessMatrix = FEMData->GetStiffnessMatrix();
	unsigned int* ColumnHeights = StiffnessMatrix->GetColumnHeights();

	for (unsigned int col = 0; col < NEQ; col++)
	{
		if (col + 1 % 10 == 0)
		{
			*this << endl;
		}

		*this << setw(8) << ColumnHeights[col];
	}

	*this << endl
		  << endl;
}

//	Print address of diagonal elements for debuging
void COutputter::PrintDiagonalAddress()
{
	*this << "*** _Debug_ *** Address of Diagonal Element" << endl;

	CDomain* FEMData = CDomain::GetInstance();

	unsigned int NEQ = FEMData->GetNEQ();
	CSkylineMatrix<double> *StiffnessMatrix = FEMData->GetStiffnessMatrix();
	unsigned int* DiagonalAddress = StiffnessMatrix->GetDiagonalAddress();

	for (unsigned int col = 0; col <= NEQ; col++)
	{
		if (col + 1 % 10 == 0)
		{
			*this << endl;
		}

		*this << setw(8) << DiagonalAddress[col];
	}

	*this << endl
		  << endl;
}

//	Print banded and full stiffness matrix for debuging
void COutputter::PrintStiffnessMatrix()
{
	*this << "*** _Debug_ *** Banded stiffness matrix" << endl;

	CDomain* FEMData = CDomain::GetInstance();

	unsigned int NEQ = FEMData->GetNEQ();
	CSkylineMatrix<double> *StiffnessMatrix = FEMData->GetStiffnessMatrix();
	unsigned int* DiagonalAddress = StiffnessMatrix->GetDiagonalAddress();

	*this << setiosflags(ios::scientific) << setprecision(5);

	for (unsigned int i = 0; i < DiagonalAddress[NEQ] - DiagonalAddress[0]; i++)
	{
		*this << setw(14) << (*StiffnessMatrix)(i);

		if ((i + 1) % 6 == 0)
		{
			*this << endl;
		}
	}

	*this << endl
		  << endl;

	*this << "*** _Debug_ *** Full stiffness matrix" << endl;

	for (int I = 1; I <= NEQ; I++)
	{
		for (int J = 1; J <= NEQ; J++)
		{
			int H = DiagonalAddress[J] - DiagonalAddress[J - 1];
			if (J - I - H >= 0)
			{
				*this << setw(14) << 0.0;
			}
			else
			{
				*this << setw(14) << (*StiffnessMatrix)(I, J);
			}
		}

		*this << endl;
	}

	*this << endl;
}

//	Print displacement vector for debuging
void COutputter::PrintDisplacement()
{
	*this << "*** _Debug_ *** Displacement vector" << endl;

	CDomain* FEMData = CDomain::GetInstance();

	unsigned int NEQ = FEMData->GetNEQ();
	double* Force = FEMData->GetForce();

	*this << setiosflags(ios::scientific) << setprecision(5);

	for (unsigned int i = 0; i < NEQ; i++)
	{
		if ((i + 1) % 6 == 0)
		{
			*this << endl;
		}

		*this << setw(14) << Force[i];
	}

	*this << endl
		  << endl;
}

#endif
