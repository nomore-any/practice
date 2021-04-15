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

//	���һ������
template <class type> void clear( type* a, unsigned int N )
{
	for (unsigned int i = 0; i < N; i++)
		a[i] = 0;
}

CDomain* CDomain::_instance = nullptr;

//	���캯������ʼ��
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

//	��������
CDomain::~CDomain()
{
	delete [] NodeList;

	delete [] EleGrpList;

	delete [] NLOAD;
	delete [] LoadCases;

	delete [] Force;
	delete StiffnessMatrix;
}

//	���_instance��һ����ָ�룬˵����û��CDomain���󣬴�ʱ�ڶ�������һ��CDomain���󣬷���ָ��ö����ָ��
CDomain* CDomain::GetInstance()
{
	if (!_instance) 
		_instance = new CDomain();
	
	return _instance;
}

//	��һ������Ĳ�����ϣ����ȡ���ļ������ڶ�������Ĳ�����ϣ����д����ļ���(outputterʵ������Ŀ�꣬һ����ֱ�Ӵ�ӡ��һ���ǰѽ��д��һ���ļ���)
bool CDomain::ReadData(string FileName, string OutFile)
{
	//��Input��ȡ�ļ�
	Input.open(FileName);
	//�ж��ļ��Ƿ���ڣ����ù�
	if (!Input) 
	{
		cerr << "*** Error *** File " << FileName << " does not exist !" << endl;
		exit(3);
	}
	//ʹ��COutputter�����GetInstance����������һ��Outputter��������������еĶ�����������д���ļ��Ľӿڣ�����ָ��
	COutputter* Output = COutputter::GetInstance(OutFile);

//	��ȡ����
	Input.getline(Title, 256);
	Output->OutputHeading();

//	��ȡ���¸���ȫ����Ϣ
	Input >> NUMNP >> NUMEG >> NLCASE >> MODEX;

//	��CDomain���ReadNodalPoints������ȡ�ڵ���Ϣ�����Ҵ�ӡ����
	if (ReadNodalPoints())
        Output->OutputNodeInfo();
    else
        return false;

//	��CDomain��ĺ�������NEQ����Outputter�����NEQ
	CalculateEquationNumber();
	Output->OutputEquationNumber();

//	��CDomain��ĺ�����ȡ�غ���Ϣ����Outputter�����
	if (ReadLoadCases())
        Output->OutputLoadInfo();
    else
        return false;

//	��CDomain��ĺ�����ȡ��Ԫ��Ϣ����Outputter�����
	if (ReadElements())
        Output->OutputElementInfo();
    else
        return false;

	return true;
}

//	��ȡ�ڵ���Ϣ
bool CDomain::ReadNodalPoints()
{

//	�ڶ�������һ������Ϊ'�ڵ�����'��CNode��������
	NodeList = new CNode[NUMNP];

//	�����еĽڵ�ѭ��
	for (unsigned int np = 0; np < NUMNP; np++)
    {
		//��һ����CNode���е�Read������ȡ���ݣ����𵽶�ȡ�����ã�Ҳ���ж������ʽ�Ƿ���ȷ������
		if (!NodeList[np].Read(Input))
			return false;
        //�жϽڵ��Ƿ��սڵ��ŵ�˳������
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

//	����������������ã�һ�Ǽ����ܷ����������ǰ�bcode��ֵ������ɶ�ȫ�ֱ��
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

//	��ȡ�غ���Ϣ
bool CDomain::ReadLoadCases()
{
//	�����غɹ�������CLoadCaseData��������
	LoadCases = new CLoadCaseData[NLCASE];	// List all load cases

//	���غɹ������ѭ��������CLoadCaseData���е�Read������ȡ�ù������غ���Ϣ
	for (unsigned int lcase = 0; lcase < NLCASE; lcase++)
    {
        unsigned int LL;
        Input >> LL;
        
		//����������˳�������غɻᱨ��
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

// ��ȡ��Ԫ��Ϣ
bool CDomain::ReadElements()
{
	//���䵥Ԫ����������
    EleGrpList = new CElementGroup[NUMEG];

//	�Ե�Ԫ����ѭ��
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
		//��ElementGroup���е�Read������ȡ��Ԫ���ݣ�ͬʱ�ж������ʽ�Ƿ���ȷ
        if (!EleGrpList[EleGrp].Read(Input))
            return false;
    
    return true;
}

//	����ȫ�ָն�����и�
void CDomain::CalculateColumnHeights()
{
#ifdef _DEBUG_
    COutputter* Output = COutputter::GetInstance();
    *Output << setw(9) << "Ele = " << setw(22) << "Location Matrix" << endl;
#endif

	//˼·�ǣ���ElementGroup��ʼѭ�����ٶ�Elementѭ������Element�Ĳ�������CSkylineMatrix��CalculateColumnHeight����������и�����
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)		//	Loop over for all element groups
    {
        CElementGroup& ElementGrp = EleGrpList[EleGrp];
        unsigned int NUME = ElementGrp.GetNUME();
        
		for (unsigned int Ele = 0; Ele < NUME; Ele++)	//	Loop over for all elements in group EleGrp
        {
            CElement& Element = ElementGrp[Ele];

            // ����һ����ÿ����Ԫ���ɷ�������(LocationMatrix���ֲ����ɶȺ�ȫ�����ɶȵı����ϵ)
            Element.GenerateLocationMatrix();
            
#ifdef _DEBUG_
            unsigned int* LocationMatrix = Element.GetLocationMatrix();
            
            *Output << setw(9) << Ele+1;
            for (int i=0; i<Element.GetND(); i++)
                *Output << setw(5) << LocationMatrix[i];
            *Output << endl;
#endif
			//�����иߣ����StiffnessMatrix��ColumnHeights_������
            StiffnessMatrix->CalculateColumnHeight(Element.GetLocationMatrix(), Element.GetND());
        }
    }
    //�����������
    StiffnessMatrix->CalculateMaximumHalfBandwidth();
    
#ifdef _DEBUG_
    *Output << endl;
	Output->PrintColumnHeights();
#endif

}

//    Ϊ���־�������������ڴ�ռ䣬���ҵ���CalculateColumnHeight������ȡ�иߣ����һЩ��ն����йص���Ϣ
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

//	��װȫ�ָն���
void CDomain::AssembleStiffnessMatrix()
{
//	˼·�ǣ��Ե�Ԫ��ѭ������ÿ����Ԫ���ڶԵ�Ԫѭ������ÿ����Ԫ����ElementStiffness������õ�Ԫ�ն��󣬴���StiffnessMatrix��Assembly�����ʵ����װ
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
		//��Ԫ��ѭ��
        CElementGroup& ElementGrp = EleGrpList[EleGrp];
        unsigned int NUME = ElementGrp.GetNUME();

		unsigned int size = ElementGrp[0].SizeOfStiffnessMatrix();
		double* Matrix = new double[size];

//		��Ԫѭ��
		for (unsigned int Ele = 0; Ele < NUME; Ele++)
        {
            CElement& Element = ElementGrp[Ele];
			//���㵥Ԫ�ն���浽Matrix������װ��ȫ�ָն���
            Element.ElementStiffness(Matrix);
            StiffnessMatrix->Assembly(Matrix, Element.GetLocationMatrix(), Element.GetND());
        }

		//��������ڴ棬�Ȳ���
		delete[] Matrix;
		Matrix = nullptr;
	}

#ifdef _DEBUG_
	COutputter* Output = COutputter::GetInstance();
	Output->PrintStiffnessMatrix();
#endif

}

//	����һ��������ţ����ɸù�����ȫ��������
bool CDomain::AssembleForce(unsigned int LoadCase)
{
	if (LoadCase > NLCASE) 
		return false;

	CLoadCaseData* LoadData = &LoadCases[LoadCase - 1];

	//��ͬ�����ǲ����ӵģ�ÿ���غɹ������ᵥ�����λ�ƣ����������Force���㣬�ټ�����һ������
    clear(Force, NEQ);

//	�Ը��غɹ������غɱ�Ž���ѭ��
	for (unsigned int lnum = 0; lnum < LoadData->nloads; lnum++)
	{
		//�ҵ���lnum���غ��������ɶȵ�ȫ�ֱ��
		//������Domain�ṩNodeList,LoadData�ṩ�ڵ��ţ��ҵ��ڵ㣬Ȼ���ýڵ��bcode�ҵ����ɶȵ�ȫ�ֱ��
		unsigned int dof = NodeList[LoadData->node[lnum] - 1].bcode[LoadData->dof[lnum] - 1];
        
        if(dof) // The DOF is activated
            Force[dof - 1] += LoadData->load[lnum];
	}

	return true;
}

