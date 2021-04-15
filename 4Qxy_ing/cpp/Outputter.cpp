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

//	��ӡʱ�䣬���ù�
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
//��COutputter�������_instance���г�ʼ������Ϊ������Բ��������κ�һ������ģ����Բ���д�ڹ��캯������ǵ�����ʼ��
COutputter* COutputter::_instance = nullptr;

//	�����ļ�������ofstream������ļ�
COutputter::COutputter(string FileName)
{
	OutputFile.open(FileName);

	if (!OutputFile)
	{
		cerr << "*** Error *** File " << FileName << " does not exist !" << endl;
		exit(3);
	}
}

//	��_instance�ǿ�ָ���ʱ������ָ��һ���½��Ķ�����COutputter���󣬲����ظ�ָ��
//  ���ϣ��������ĺ�������Outputter���������ְ취��һ�־���ʹ�ø�ָ�룬��һ�־���ֱ�Ӵ���outputter���󣬻����϶��Ǻ���
COutputter* COutputter::GetInstance(string FileName)
{
	if (!_instance)
		_instance = new COutputter(FileName);
    
	return _instance;
}

//	���ù�
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

//	��ӡ�ڵ���Ϣ
void COutputter::OutputNodeInfo()
{
	//��CDomain��ȡָ��CDomain�����ָ��
	CDomain* FEMData = CDomain::GetInstance();
	//��CDomain��ȡCNode���͵�����NodeList�����ڵ�����
	CNode* NodeList = FEMData->GetNodeList();

	//���ù�
	*this << "C O N T R O L   I N F O R M A T I O N" << endl
		  << endl;
	*this << setiosflags(ios::scientific) << setprecision(5);

	//����FEMDataָ���ȡCDomain�����е�'�ڵ�����'��'��Ԫ�����'��'�غɹ�����'��'���ģʽ'
	unsigned int NUMNP = FEMData->GetNUMNP();
	unsigned int NUMEG = FEMData->GetNUMEG();
	unsigned int NLCASE = FEMData->GetNLCASE();
	unsigned int MODEX = FEMData->GetMODEX();

	//��ӡ'�ڵ�����'��'��Ԫ�����'��'�غɹ�����'��'���ģʽ'
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

	//����CNode���е�Write�����Ѹ������ݴ�ӡ�����Լ�д��Ŀ���ļ���
	for (unsigned int np = 0; np < NUMNP; np++)
		NodeList[np].Write(*this);

	*this << endl;
}

//	�Խڵ��Ž���ѭ�������ÿ���ڵ��3�����ɶȶ�Ӧ��ȫ�ֱ��(EquationNumber)
void COutputter::OutputEquationNumber()
{
	//��CDomain��ȡָ��CDomain�����ָ��
	CDomain* FEMData = CDomain::GetInstance();
	//ʹ��FEMDataָ���ȡ'�ڵ�����'
	unsigned int NUMNP = FEMData->GetNUMNP();
	//ʹ��FEMDataָ���ȡCDomain������CNode���͵�����NodeList
	CNode* NodeList = FEMData->GetNodeList();

	*this << " EQUATION NUMBERS" << endl
		  << endl;
	*this << "   NODE NUMBER   DEGREES OF FREEDOM" << endl;
	*this << "        N           X    Y    Z" << endl;
	//�Խڵ�ѭ��������CNode���е�WriteEquationNo����ýڵ��������ɶȵ�ȫ�ֱ��
	for (unsigned int np = 0; np < NUMNP; np++) // Loop over for all node
		NodeList[np].WriteEquationNo(*this);

	*this << endl;
}

//	�����Ԫ��Ϣ
void COutputter::OutputElementInfo()
{
	//��ȡָ��CDomain�����ָ��
	CDomain* FEMData = CDomain::GetInstance();
	//��CDomain�����л�ȡ��Ԫ����
	unsigned int NUMEG = FEMData->GetNUMEG();

	*this << " E L E M E N T   G R O U P   D A T A" << endl
		  << endl
		  << endl;
	//�Ե�Ԫ���Ž���ѭ��
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
		*this << " E L E M E N T   D E F I N I T I O N" << endl
			  << endl;
		//Ϊ�˻�ȡ�õ�Ԫ��ĵ�Ԫ���ͣ��ȴӵ�Ԫ������õ���Ԫ�飬���õ�Ԫ�����е�GetElementType()������ȡ��Ԫ����
		//ͬ��Ҳ�����ַ�����ȡ��Ԫ���еĵ�Ԫ����NUME
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
				//�������������档�����ElementGroup�����е�Ԫ����Ϣ
				OutputBarElements(EleGrp);
				break;
			case ElementTypes::Q4: // Q4 element
				//�������������档�����ElementGroup�����е�Ԫ����Ϣ
				OutputQ4Elements(EleGrp);
				break;
		    default:
		        *this << ElementType << " has not been implemented yet." << endl;
		        break;
		}
	}
}
//    ���뵥Ԫ����(�����Ǹ˵�Ԫ�ĵ�Ԫ��)������õ�Ԫ�����е�Ԫ����Ϣ
void COutputter::OutputBarElements(unsigned int EleGrp)
{
	//��ȡָ��CDomain�����ָ��
	CDomain* FEMData = CDomain::GetInstance();
	//��ȡ�õ�Ԫ������Լ���Ԫ��Ĳ�������
	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();
	//���������Ϣ
	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		  << endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
	*this << " AND CROSS-SECTIONAL  CONSTANTS  . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
		  << endl
		  << endl;

	*this << "  SET       YOUNG'S     CROSS-SECTIONAL" << endl
		  << " NUMBER     MODULUS          AREA" << endl
		  << "               E              A" << endl;
	//û����@
	*this << setiosflags(ios::scientific) << setprecision(5);

	for (unsigned int mset = 0; mset < NUMMAT; mset++)
    {
        *this << setw(5) << mset+1;
	//������ϲ����ķ�����ʹ��ElementGroup���е�GetMaterial�����ҵ���Ԫ����ĳ�ֲ��ϵĶ�������CBarMaterial���е�Write�����������
		ElementGroup.GetMaterial(mset).Write(*this);
    }

	*this << endl << endl
		  << " E L E M E N T   I N F O R M A T I O N" << endl;
    
	*this << " ELEMENT     NODE     NODE       MATERIAL" << endl
		  << " NUMBER-N      I        J       SET NUMBER" << endl;
	//������е�Ԫ��Ϣ��ԭ��ͬ��
	unsigned int NUME = ElementGroup.GetNUME();

	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < NUME; Ele++)
    {
        *this << setw(5) << Ele+1;
		ElementGroup[Ele].Write(*this);
    }

	*this << endl;
}

//    ���뵥Ԫ����(������4Q��Ԫ�ĵ�Ԫ��)������õ�Ԫ�����е�Ԫ����Ϣ
void COutputter::OutputQ4Elements(unsigned int EleGrp)
{
	//��ȡָ��CDomain�����ָ��
	CDomain* FEMData = CDomain::GetInstance();
	//��ȡ�õ�Ԫ������Լ���Ԫ��Ĳ�������
	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();
	//���������Ϣ
	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		  << endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
	*this << " AND CROSS-SECTIONAL  CONSTANTS  . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
		  << endl
		  << endl;

	*this << "  SET       YOUNG'S     CROSS-SECTIONAL" << endl
		  << " NUMBER     MODULUS          ���ɱ�" << endl
		  << "               E              Posion" << endl;
	//û����@
	*this << setiosflags(ios::scientific) << setprecision(5);

	for (unsigned int mset = 0; mset < NUMMAT; mset++)
    {
        *this << setw(5) << mset+1;
	//������ϲ����ķ�����ʹ��ElementGroup���е�GetMaterial�����ҵ���Ԫ����ĳ�ֲ��ϵĶ�������CQ4Material���е�Write�����������
		ElementGroup.GetMaterial(mset).Write(*this);
    }

	*this << endl << endl
		  << " E L E M E N T   I N F O R M A T I O N" << endl;
    
	*this << " ELEMENT     NODE     NODE     NODE     NODE     MATERIAL" << endl
		  << " NUMBER-N     I        J        K        L      SET NUMBER" << endl;
	//������е�Ԫ��Ϣ��ԭ��ͬ��
	unsigned int NUME = ElementGroup.GetNUME();

	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < NUME; Ele++)
    {
        *this << setw(5) << Ele+1;
		ElementGroup[Ele].Write(*this);
    }

	*this << endl;
}

//	��ӡ�غ���Ϣ
void COutputter::OutputLoadInfo()
{
	//��ȡָ��CDomain�����ָ��
	CDomain* FEMData = CDomain::GetInstance();
	//�Թ������ѭ��
	for (unsigned int lcase = 1; lcase <= FEMData->GetNLCASE(); lcase++)
	{
		//��ȡָ��ĳ��CLoadCaseData�����ָ��
		CLoadCaseData* LoadData = &FEMData->GetLoadCases()[lcase - 1];

		//������ֶ���
		*this << setiosflags(ios::scientific);
		*this << " L O A D   C A S E   D A T A" << endl
			  << endl;

		*this << "     LOAD CASE NUMBER . . . . . . . =" << setw(6) << lcase << endl;
		*this << "     NUMBER OF CONCENTRATED LOADS . =" << setw(6) << LoadData->nloads << endl
			  << endl;
		*this << "    NODE       DIRECTION      LOAD" << endl
			  << "   NUMBER                   MAGNITUDE" << endl;
		//ʹ��CLoadCaseData���µ�Write�������
		LoadData->Write(*this);

		*this << endl;
	}
}

//	���λ����Ϣ
void COutputter::OutputNodalDisplacement()
{
	//��CDomain���ȡָ�롢�ڵ��б�λ������
	CDomain* FEMData = CDomain::GetInstance();
	CNode* NodeList = FEMData->GetNodeList();
	double* Displacement = FEMData->GetDisplacement();

	*this << setiosflags(ios::scientific);

	*this << " D I S P L A C E M E N T S" << endl
		  << endl;
	*this << "  NODE           X-DISPLACEMENT    Y-DISPLACEMENT    Z-DISPLACEMENT" << endl;
	//���ÿ���ڵ��λ����Ϣ��ʹ��CNode���µ�WriteNodalDisplacement����
	for (unsigned int np = 0; np < FEMData->GetNUMNP(); np++)
		NodeList[np].WriteNodalDisplacement(*this, Displacement);

	*this << endl;
}

//	�����ԪӦ��
void COutputter::OutputElementStress()
{
	//��ȡָ��
	CDomain* FEMData = CDomain::GetInstance();
	//��ȡȫ��λ������
	double* Displacement = FEMData->GetDisplacement();
	//��ȡ��Ԫ�����
	unsigned int NUMEG = FEMData->GetNUMEG();
    //�Ե�Ԫ��ѭ��
	for (unsigned int EleGrpIndex = 0; EleGrpIndex < NUMEG; EleGrpIndex++)
	{
		*this << " S T R E S S  C A L C U L A T I O N S  F O R  E L E M E N T  G R O U P" << setw(5)
			  << EleGrpIndex + 1 << endl
			  << endl;
		//��CDomain���ȡ��Ԫ��
		CElementGroup& EleGrp = FEMData->GetEleGrpList()[EleGrpIndex];
		//�ӵ�Ԫ���ȡ��Ԫ��������Ԫ����
		unsigned int NUME = EleGrp.GetNUME();
		ElementTypes ElementType = EleGrp.GetElementType();

		switch (ElementType)
		{
			case ElementTypes::Bar: // Bar element
				*this << "  ELEMENT             FORCE            STRESS" << endl
					<< "  NUMBER" << endl;

				double stress;
				//�Ե�Ԫѭ�������ÿ����Ԫ��Ӧ��
				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					CElement& Element = EleGrp[Ele];
					Element.ElementStress(&stress, Displacement);
		//����GetElementMaterial()��д�ڸ�������ģ�����ָ��CMaterial�࣬��Ҫ����ǿ������ת�����ѽ����Ϊָ��CBarMaterial��
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
				//�Ե�Ԫѭ�������ÿ����Ԫ��Ӧ��
				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					CElement& Element = EleGrp[Ele];
					Element.ElementStress(&stress, Displacement);
		//����GetElementMaterial()��д�ڸ�������ģ�����ָ��CMaterial�࣬��Ҫ����ǿ������ת�����ѽ����Ϊָ��CBarMaterial��
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

//	���һЩ��ն����йصĶ���
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
