#pragma once

#include "Element.h"

using namespace std;

//! ��Element��̳е�Bar��
class CQ4 : public CElement
{
public:

//!	���캯��
	CQ4();

//!	��������
	~CQ4();

//!	��ȡQ4��Ԫ������
	virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

//!	��ӡQ4��Ԫ������
	virtual void Write(COutputter& output);

//!	���㵥Ԫ�ն��󣬴浽�����������
	virtual void ElementStiffness(double* Matrix);

//!	���㵥ԪӦ��
	virtual void ElementStress(double* stress, double* Displacement);
};
