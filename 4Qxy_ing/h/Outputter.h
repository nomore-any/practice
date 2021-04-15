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
#include <iostream>
#include <iomanip>

using namespace std;

//! Outputer class is used to output results
class COutputter
{
private:

//!	����һ��ofstream��Ķ���
	ofstream OutputFile;

//!	����ָ��COutputter���ָ�룬���ָ���Ǿ�̬�ģ���˲��Ƕ�������ԣ����������������
	static COutputter* _instance;

//! COutputter�౻����ʱ��һ���ļ��������Ҫһ���ļ���������
    COutputter(string FileName);

public:

//!	����ofstream������ָ��
	inline ofstream* GetOutputFile() { return &OutputFile; }

//!	������ĳ��COutputter����������ļ���
	static COutputter* GetInstance(string FileName = " ");

//!	Output current time and date
	void PrintTime(const struct tm * ptm, COutputter& output);

//!	���ù�
	void OutputHeading();

//!	��ӡ�ڵ���Ϣ
	void OutputNodeInfo();

//!	��ӡ��������Ϣ
	void OutputEquationNumber();

//!	��ӡ��Ԫ��Ϣ
	void OutputElementInfo();

//!	��ӡ��Ԫ����Ϣ
	void OutputBarElements(unsigned int EleGrp);
	void OutputQ4Elements(unsigned int EleGrp);

//!	��ӡ�غ���Ϣ
	void OutputLoadInfo(); 

//!	��ӡλ����Ϣ
	void OutputNodalDisplacement();

//!	��ӡӦ����Ϣ
	void OutputElementStress();

//!	��ӡ����ϵͳ����Ϣ
	void OutputTotalSystemData();

//! ������<<������������������͵����ݣ�ͬʱ��ӡ���ݺ����ļ���д�����룻����COutputter�����Ա��������
	template <typename T>
	COutputter& operator<<(const T& item) 
	{
		std::cout << item;
		OutputFile << item;
		return *this;
	}
//! ����������������أ�����û�п���@
	typedef std::basic_ostream<char, std::char_traits<char> > CharOstream;
	COutputter& operator<<(CharOstream& (*op)(CharOstream&)) 
	{
		op(std::cout);
		op(OutputFile);
		return *this;
	}

#ifdef _DEBUG_

//!	Print banded and full stiffness matrix for debuging
	void PrintStiffnessMatrix();

//!	Print address of diagonal elements for debuging
	void PrintDiagonalAddress();

//!	Print column heights for debuging
	void PrintColumnHeights();

//!	Print displacement vector for debuging
	void PrintDisplacement();

#endif

};
