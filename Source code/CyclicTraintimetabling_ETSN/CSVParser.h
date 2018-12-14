#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <sstream>
using namespace std;
using std::string;
using std::ifstream;
using std::vector;
using std::map;
using std::istringstream;
template <typename T>

#pragma warning(disable: 4244)  // stop warning: "conversion from 'int' to 'float', possible loss of data"

string NumberToString(T Number)
{
	ostringstream ss;
	ss << Number;
	return ss.str();
}


template <typename T>
T StringToNumber(const string &Text)
{
	istringstream ss(Text);
	T result;
	return ss >> result ? result : 0;
}
class CCSVParser
{
public:
	char Delimiter;
	bool IsFirstLineHeader;
	ifstream inFile;
	vector<string> LineFieldsValue;
	vector<string> Headers;
	map<string, int> FieldsIndices;

	vector<int> LineIntegerVector;

public:
	void  ConvertLineStringValueToIntegers()
	{
		LineIntegerVector.clear();
		for (unsigned i = 0; i < LineFieldsValue.size(); i++)
		{
			std::string si = LineFieldsValue[i];
			int value = atoi(si.c_str());

			if (value >= 1)
				LineIntegerVector.push_back(value);

		}
	}
	vector<string> GetHeaderVector()
	{
		return Headers;
	}

	int m_EmptyLineCount;
	bool m_bDataHubSingleCSVFile;
	string m_DataHubSectionName;
	bool m_bLastSectionRead;

	bool m_bSkipFirstLine;  // for DataHub CSV files

	CCSVParser::CCSVParser(void)
	{
		Delimiter = ',';
		IsFirstLineHeader = true;
		m_bSkipFirstLine = false;
		m_bDataHubSingleCSVFile = false;
		m_bLastSectionRead = false;
		m_EmptyLineCount++;
	}

	CCSVParser::~CCSVParser(void)
	{
		if (inFile.is_open()) inFile.close();
	}


	bool CCSVParser::OpenCSVFile(string fileName, bool bIsFirstLineHeader)
	{
		inFile.clear();
		inFile.open(fileName.c_str());

		IsFirstLineHeader = bIsFirstLineHeader;
		if (inFile.is_open())
		{
			if (m_bSkipFirstLine)
			{
				string s;
				std::getline(inFile, s);
			}
			if (IsFirstLineHeader)
			{
				string s;
				std::getline(inFile, s);

				if (s.length() == 0)
					return true;

				vector<string> FieldNames = ParseLine(s);

				for (size_t i = 0; i<FieldNames.size(); i++)
				{
					string tmp_str = FieldNames.at(i);
					size_t start = tmp_str.find_first_not_of(" ");

					string name;
					if (start == string::npos)
					{
						name = "";
					}
					else
					{
						name = tmp_str.substr(start);
						TRACE("%s,", name.c_str());
					}
					Headers.push_back(name);
					FieldsIndices[name] = (int)i;
				}
			}

			return true;
		}
		else
		{
			return false;
		}
	}

	void CCSVParser::CloseCSVFile(void)
	{
		inFile.close();
	}



	bool CCSVParser::ReadRecord()
	{
		LineFieldsValue.clear();

		if (inFile.is_open())
		{
			string s;
			std::getline(inFile, s);
			if (s.length() > 0)
			{
				LineFieldsValue = ParseLine(s);
				return true;
			}
			else
			{
				return false;
			}
		}
		else
		{
			return false;
		}
	}

	vector<string> CCSVParser::ParseLine(string line)
	{
		vector<string> SeperatedStrings;
		string subStr;

		if (line.length() == 0)
			return SeperatedStrings;

		istringstream ss(line);


		if (line.find_first_of('"') == string::npos)
		{

			while (std::getline(ss, subStr, Delimiter))
			{
				SeperatedStrings.push_back(subStr);
			}

			if (line.at(line.length() - 1) == ',')
			{
				SeperatedStrings.push_back("");
			}
		}
		else
		{
			while (line.length() > 0)
			{
				size_t n1 = line.find_first_of(',');
				size_t n2 = line.find_first_of('"');

				if (n1 == string::npos && n2 == string::npos) //last field without double quotes
				{
					subStr = line;
					SeperatedStrings.push_back(subStr);
					break;
				}

				if (n1 == string::npos && n2 != string::npos) //last field with double quotes
				{
					size_t n3 = line.find_first_of('"', n2 + 1); // second double quote

					//extract content from double quotes
					subStr = line.substr(n2 + 1, n3 - n2 - 1);
					SeperatedStrings.push_back(subStr);

					break;
				}

				if (n1 != string::npos && (n1 < n2 || n2 == string::npos))
				{
					subStr = line.substr(0, n1);
					SeperatedStrings.push_back(subStr);
					if (n1 < line.length() - 1)
					{
						line = line.substr(n1 + 1);
					}
					else // comma is the last char in the line string, push an empty string to the back of vector
					{
						SeperatedStrings.push_back("");
						break;
					}
				}

				if (n1 != string::npos && n2 != string::npos && n2 < n1)
				{
					size_t n3 = line.find_first_of('"', n2 + 1); // second double quote
					subStr = line.substr(n2 + 1, n3 - n2 - 1);
					SeperatedStrings.push_back(subStr);
					size_t idx = line.find_first_of(',', n3 + 1);

					if (idx != string::npos)
					{
						line = line.substr(idx + 1);
					}
					else
					{
						break;
					}
				}
			}

		}

		return SeperatedStrings;
	}
	
	template <class T> bool GetValueByFieldName(string field_name, T& value, bool NonnegativeFlag = true)
	{
		if (FieldsIndices.find(field_name) == FieldsIndices.end())
		{
			return false;
		}
		else
		{
			if (LineFieldsValue.size() == 0)
			{
				return false;
			}

			int size = (int)(LineFieldsValue.size());
			if (FieldsIndices[field_name] >= size)
			{
				return false;
			}

			string str_value = LineFieldsValue[FieldsIndices[field_name]];

			if (str_value.length() <= 0)
			{
				return false;
			}

			istringstream ss(str_value);

			T converted_value;
			ss >> converted_value;

			if (/*!ss.eof() || */ ss.fail())
			{
				return false;
			}

			if (NonnegativeFlag && converted_value<0)
				converted_value = 0;

			value = converted_value;
			return true;
		}
	}

	bool GetValueByFieldName(string field_name, string& value)
	{
		if (FieldsIndices.find(field_name) == FieldsIndices.end())
		{
			return false;
		}
		else
		{
			if (LineFieldsValue.size() == 0)
			{
				return false;
			}

			unsigned int index = FieldsIndices[field_name];
			if (index >= LineFieldsValue.size())
			{
				return false;
			}
			string str_value = LineFieldsValue[index];

			if (str_value.length() <= 0)
			{
				return false;
			}

			value = str_value;
			return true;
		}
	}

	template <class T> bool GetValueBySectionKeyFieldName(string file_name, string section_name, string key_name, string field_name, T& value)
	{
		OpenCSVFile(file_name, true);

		while (ReadRecord())
		{
			if (LineFieldsValue[0] != section_name || LineFieldsValue[1] != key_name)
				continue;

			if (FieldsIndices.find(field_name) == FieldsIndices.end())
			{
				CloseCSVFile();
				return false;
			}
			else
			{
				if (LineFieldsValue.size() == 0)
				{
					CloseCSVFile();
					return false;
				}

				int size = (int)(LineFieldsValue.size());
				if (FieldsIndices[field_name] >= size)
				{
					CloseCSVFile();
					return false;
				}

				string str_value = LineFieldsValue[FieldsIndices[field_name]];

				if (str_value.length() <= 0)
				{
					CloseCSVFile();
					return false;
				}

				istringstream ss(str_value);

				T converted_value;
				ss >> converted_value;

				if (/*!ss.eof() || */ ss.fail())
				{

					CloseCSVFile();
					return false;
				}

				value = converted_value;
				CloseCSVFile();
				return true;
			}
		}

		CloseCSVFile();

		return false;
	}

};

