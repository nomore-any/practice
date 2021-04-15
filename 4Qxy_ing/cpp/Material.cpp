/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Material.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

//	读取材料参数
bool CBarMaterial::Read(ifstream& Input)
{
	Input >> nset;	//读数材料编号(每个单元组都有自己从1开始的材料编号)

	Input >> E >> Area;	// 读取E和A

	return true;
}

//	输出材料参数
void CBarMaterial::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << Area << endl;
}


//	读取材料参数
bool CQ4Material::Read(ifstream& Input)
{
	Input >> nset;	//读数材料编号(每个单元组都有自己从1开始的材料编号)

	Input >> E >> Nu;	// 读取E和A和泊松比

	return true;
}

//	输出材料参数
void CQ4Material::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << Nu << endl;
}
