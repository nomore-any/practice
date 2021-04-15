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

//!	创建一个ofstream类的对象
	ofstream OutputFile;

//!	创建指向COutputter类的指针，这个指针是静态的，因此不是对象的属性，而是整个类的属性
	static COutputter* _instance;

//! COutputter类被创建时打开一个文件，因此需要一个文件名的输入
    COutputter(string FileName);

public:

//!	传回ofstream类对象的指针
	inline ofstream* GetOutputFile() { return &OutputFile; }

//!	传回与某个COutputter对象关联的文件名
	static COutputter* GetInstance(string FileName = " ");

//!	Output current time and date
	void PrintTime(const struct tm * ptm, COutputter& output);

//!	不用管
	void OutputHeading();

//!	打印节点信息
	void OutputNodeInfo();

//!	打印方程数信息
	void OutputEquationNumber();

//!	打印单元信息
	void OutputElementInfo();

//!	打印单元组信息
	void OutputBarElements(unsigned int EleGrp);
	void OutputQ4Elements(unsigned int EleGrp);

//!	打印载荷信息
	void OutputLoadInfo(); 

//!	打印位移信息
	void OutputNodalDisplacement();

//!	打印应力信息
	void OutputElementStress();

//!	打印整个系统的信息
	void OutputTotalSystemData();

//! 重载了<<运算符，输入任意类型的数据，同时打印数据和在文件中写入输入；返回COutputter对象，以便连续输出
	template <typename T>
	COutputter& operator<<(const T& item) 
	{
		std::cout << item;
		OutputFile << item;
		return *this;
	}
//! 发生了运算符的重载，但是没有看懂@
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
