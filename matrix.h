//matrix.h
#ifndef MAT_H_H
#define MAT_H_H
#include<iostream>
#include<cstdlib>
#include<stdexcept>
#include<exception>

class Matrix
{
private:
	int row = 1;
	int col = 1;
	double* item;
public:
	Matrix() : row(0), col(0), item(nullptr) {}
	Matrix(const int x ,const int y,const double z= 0.0) :row(x),col(y)//���캯��
	{
		item = new double[x*y];//���䶯̬�ڴ�
		for (int i = 0;i < x*y;i++)
		{
			item[i] = z;
		}
	}
	//����ʹ�ö�̬�ڴ���࣬������ʾ�Ķ����乹�캯�����������������ƹ��캯�������ظ�ֵ���������������ͷ��ڴ���ٴη��ʣ���ɷ����޶�����ڴ棬���ͷŶ��ͬһ���ڴ���ɳ������������ʾ����ĸ��ƹ��캯�������ظ�ֵ������У����������ȸ��ƣ�ҪΪ�����·����µ��ڴ��ַ���ٰ����ݸ��ƹ�����
	//����ô���ĺ���Ǳ���������Ĭ�ϵĸ��ƹ��캯����ֵ�����������������ָ���Ա����ָ��ͬһ���ڴ��ַ����һ���������ʱ���������ͷŸ��ڴ棬���������Ǹ��������������ʸó�Ա���������������
	~Matrix()//��������
	{
		delete[] item;//�ͷ��ڴ�
	}
	Matrix(const Matrix& mt)//���ƹ��캯��
	{
		item = new double[mt.row * mt.col];
		row = mt.row;
		col = mt.col;
		for (int i = 0;i < mt.row*mt.col ;i++)
		{
			item[i] = mt.item[i];
		}
	}
	Matrix& operator=(const Matrix& mt)//���ظ�ֵ�����
	{
		if (this == &mt)
			return *this;
		row = mt.row;
		col = mt.col;
		delete[] item;
		item = new double[mt.row * mt.col];
		for (int i = 0;i < mt.row * mt.col ;i++)
		{
			item[i] = mt.item[i];
		}
		return *this;
	}
	bool operator==(const Matrix& mt)
	{
		if (row != mt.row || col != mt.col)
			return false;
		for (size_t i = 0; i < row; i++)
		{
			for (size_t j = 0; j < col; j++)
			{
				if (item[i * col + j] != mt.item[i * col + j])
					return false;
			}
		}
		return true;
	}
	bool operator!=(const Matrix& mt)
	{
		return !(*this == mt);
	}
	Matrix operator+(const Matrix& mt) const
	{
		Matrix R(row, col, 0.0);
		try
		{
			if (row != mt.row || col != mt.col)
				throw std::length_error("+ must be operated between two matrices with the same size");
			for (int i = 0;i < row;i++)
			{
				for (int j = 0;j < col;j++)
				{
					R.item[i * col + j] = item[i * col + j] + mt.item[i * col + j];
				}
			}
			return R;
		}
		catch (std::length_error& e)
		{
			std::cout << e.what() << std::endl;
			exit(EXIT_FAILURE);
		}
	}
	Matrix operator-(const Matrix& mt) const
	{
		Matrix R(row, col, 0.0);
		try
		{
			if (row != mt.row || col != mt.col)
				throw std::length_error("- must be operated between two matrices with the same size");
			for (int i = 0;i < row;i++)
			{
				for (int j = 0;j < col;j++)
				{
					R.item[i * col + j] = item[i * col + j] - mt.item[i * col + j];
				}
			}
			return R;
		}
		catch (std::length_error& e)
		{
			std::cout << e.what() << std::endl;
			exit(EXIT_FAILURE);
		}
	}
	Matrix operator*(const Matrix& mt) const
	{
		Matrix R(row, mt.col, 0.0);
		if (col != mt.row)
			throw std::out_of_range("left matrix's column must be equal to right matrix's row");
		for (int i = 0;i < row;i++)
		{
			for (int j = 0;j < mt.col;j++)
			{
				for (int k = 0;k < col;k++)
				{
					R.item[i * col + j] += item[i * col + k] * mt.item[k * col + j];
				}
			}
		}
		return R;
	}
	Matrix operator*(const double x) const//�������أ���*�Ҳ���doubleʱ���ã������Ԫ��ͬ�˸�ֵ
	{
		Matrix R(row, col);
		for (int i = 0;i < row;i++)
		{
			for (int j = 0;j < col;j++)
			{
				R.item[i * col + j] = item[i * col + j] * x;
			}
		}
		return R;
	}
	friend Matrix operator*(const double x,const Matrix& mt)//��Ԫ������������*doubleʱ����
	{
		Matrix R(mt.row, mt.col);
		for (int i = 0;i < mt.row;i++)
		{
			for (int j = 0;j < mt.col;j++)
			{
				R.item[i * mt.col + j] = mt.item[i * mt.col + j] * x;
			}
		}
		return R;
	}
	Matrix operator/(const double x) const //������/doubleʱ����
	{
		Matrix R(row, col);
		for (int i = 0;i < row;i++)
		{
			for (int j = 0;j < col;j++)
			{
				R.item[i * col + j] = item[i * col + j] / x;
			}
		}
		return R;
	}
	const double& operator()(const int x,const int y) const
	{
		try 
		{
			if (x < 0 || x >= row || y < 0 || y >= col)
			{
				throw std::out_of_range("Index out of bounds");
			}
			return item[x * col + y];
		}
		catch (std::out_of_range& e)
		{
			std::cout << e.what() << std::endl;
			exit(EXIT_FAILURE);
		}
	}
	double& operator()(const int x, const int y)
	{
		try
		{
			if (x < 0 || x >= row || y < 0 || y >= col)
			{
				throw std::out_of_range("Index out of bounds");
			}
			return item[x * col + y];
		}
		catch (std::out_of_range& e)
		{
			std::cout << e.what() << std::endl;
			exit(EXIT_FAILURE);
		}
	}
	 
	const int& get_xsize()const { return row; }
	const int& get_ysize()const { return col; }
	void disp()const
	{
		for (int i = 0;i < row;i++)
		{
			for (int j = 0;j < col;j++)
			{
				std::cout << item[i * col + j] << "     ";
			}
			std::cout << "\n";
		}
		std::cout << "\n";
	}

	template <size_t Rows, size_t Cols>
	// set���������ڽ��ն�ά������Ϊ�������и�ֵ������ C++ ��ģ�����ԣ��ñ��������ݴ����ʵ�ʶ�ά�����������Զ��Ƶ������ά����Ϣ
	void set(const double (&a)[Rows][Cols])
	{
		delete[] item;
		row = Rows;
		col = Cols;
		item = new double[row * col];
		for (size_t i = 0; i < row; i++)
		{
			for (size_t j = 0; j < col; j++)
			{
				item[i * col + j] = a[i][j];
			}
		}
	}

	// ����set���������ڽ��ճ�ʼ���б���Ϊ�������и�ֵ
	 // ʹ�ó�ʼ���б�ֵ����һ�ν���gpt��̫�鷳��
	void set(std::initializer_list<std::initializer_list<double>> initList)
	{
		// ��ȡ��ʼ���б���к���
		row = initList.size();
		col = initList.begin()->size();

		// ȷ�������С�ͳ�ʼ���б�ƥ��
		item = new double[row * col];
		int i = 0;
		for (auto& rowList : initList)
		{
			int j = 0;
			for (auto& val : rowList)
			{
				item[i * col + j] = val;
				j++;
			}
			i++;
		}
	}
	const Matrix T() const
	{
		//���ؾ����ת�ã����ı������
		Matrix R(col, row, 0.0);
		for (int i = 0;i < col;i++)
		{
			for (int j = 0;j < row;j++)
			{
				R.item[i * row + j] = item[j * col + i];
			}
		}
		return R;
	}
	void swap_row(int x, int y)
	{
		//�����������ĵ�x�к͵�y�У�x��y��0��ʼ����
		double temp;
		for (int j = 0;j < col;j++)
		{
			temp = item[x * col + j];
			item[x * col + j] = item[y * col + j];
			item[y * col + j] = temp;
		}
	}
	Matrix standard_foramt() const
	{
		//���ؾ���任Ϊ����������ʽ����������
		Matrix R = *this;
		bool change;
		for (int j = 0;j < col;j++)
		{
			if (R.item[j * col + j] == 0)
			{
				change = false;
				for (int k = 1;k < row - j;k++)
				{
					if (R.item[(j + k) * col + j] != 0)
					{
						R.swap_row(j, j + k);//���°����ң�ֱ��ȷ��(j,j)��Ϊ0
						change = true;
						break;
					}
				}
			}
			else
				change = true;
			if (!change)//��һ�ж��Ҳ�����Ϊ0��Ԫ�أ�ֱ�Ӵ�����һ��
				continue;
			for (int i = j + 1;i < row;i++)
			{
				//�����·�ÿһ��
				double coe = R.item[i * col + j] / R.item[j * col + j];
				for (int k = 0;k < col;k++)
				{
					R.item[i * col + k] -= coe * R.item[j * col+ k];
				}
			}
		}
		return R;
	}
	Matrix accompany() const
	{
		//���ض������İ�����󣬶�������
		try
		{
			if (row != col)
				throw std::out_of_range("Only square matrix has accompany matrix");
			if (row == 1)
			{
				Matrix R(1, 1, 1.0);
				return R;
			}
			Matrix R(row, col), Temp(row - 1, col - 1);
			for (int i = 0;i < row;i++)
			{
				for (int j = 0;j < col;j++)
				{
					//����Temp��������Ϊthis��ȥ��i�е�j��
					for (int m = 0;m < row;m++)
					{
						for (int n = 0;n < col;n++)
						{
							if (m == i || n == j)
								continue;
							//����this�ĵ�i�е�j�е�Temp�ĵ�m>i?(m-1):m�е�n>j?(n-1):n��
							Temp.item[(m > i ? (m - 1) : m) * (col - 1) + (n > j ? (n - 1) : n)] = item[m * col + n];
						}
					}
					R.item[i * col + j] = ((i + j) % 2) == 0 ? det(Temp) : -det(Temp);
				}
			}
			return R.T();
		}
		catch (std::out_of_range& e)
		{
			std::cout << e.what() << std::endl;
			exit(EXIT_FAILURE);
		}
	}
	friend double det(const Matrix& mt)
	{
		//���ض��������ʽֵ��������Ԫ����
		if (mt.row != mt.col)
			throw"The matrix's row must be equal to its column";
		Matrix R = mt;

		bool change ,row_change = true;
		for (int j = 0;j < mt.col;j++)
		{
			if (R.item[j * mt.col + j] == 0)
			{
				change = false;
				for (int k = 1;k < mt.row - j;k++)
				{
					if (R.item[(j + k) * mt.col + j] != 0)
					{
						R.swap_row(j, j + k);//���°����ң�ֱ��ȷ��(j,j)��Ϊ0
						change = true;
						row_change = !row_change;
						break;
					}
				}
			}
			else
				change = true;
			if (!change)//��һ�ж��Ҳ�����Ϊ0��Ԫ�أ�ֱ�Ӵ�����һ��
				continue;
			for (int i = j + 1;i < mt.row;i++)
			{
				//�����·�ÿһ��
				double coe = R.item[i * mt.col + j] / R.item[j * mt.col + j];
				for (int k = 0;k < mt.col;k++)
				{
					R.item[i * mt.col + k] -= coe * R.item[j * mt.col + k];
				}
			}
		}
		//�����R���Ƕ���������������ʽ����standard_foramt()������ͬ����������Ҫ�����ȡ�������в����Ĵ�����ÿ��������һ�Σ�����ʽֵȡһ�θ���
		double results = 1.0;
		for (int i = 0;i < mt.row;i++)
		{
			results *= R.item[i * mt.col + i];//�Խ�Ԫ���۳˵õ�����ʽֵ
		}
		//�ڻ������ǹ����У�ÿ����һ�ν������еĲ���������ʽҪȡһ�θ���
		//���������������ĺ����ǣ���row_changeΪ��ʱ����results����row_changeΪ��ʱ����-results
		return row_change ? results : -results;
	}
	friend Matrix inv(const Matrix& mt)
	{
		//��������󣬲�����Ԫ������ʽ
		if (mt.row != mt.col)
			throw std::range_error("Only square matrix has inverse matrix");
		//�������ڰ�������������ʽ������� / �Ѿ����أ�����ֱ��ʹ��
		return mt.accompany() / det(mt);
	}

};

#endif