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
	Matrix(const int x ,const int y,const double z= 0.0) :row(x),col(y)//构造函数
	{
		item = new double[x*y];//分配动态内存
		for (int i = 0;i < x*y;i++)
		{
			item[i] = z;
		}
	}
	//对于使用动态内存的类，必须显示的定义其构造函数、析构函数、复制构造函数、重载赋值运算符，避免程序释放内存后再次访问，造成访问无定义的内存，或释放多次同一块内存造成程序崩溃。在显示定义的复制构造函数、重载赋值运算符中，必须进行深度复制，要为其重新分配新的内存地址，再把内容复制过来。
	//不这么做的后果是编译器按照默认的复制构造函数或赋值运算符，将多个对象的指针成员变量指向同一块内存地址，当一个对象过期时程序主动释放该内存，造成另外的那个对象不能正常访问该成员变量，程序崩溃。
	~Matrix()//析构函数
	{
		delete[] item;//释放内存
	}
	Matrix(const Matrix& mt)//复制构造函数
	{
		item = new double[mt.row * mt.col];
		row = mt.row;
		col = mt.col;
		for (int i = 0;i < mt.row*mt.col ;i++)
		{
			item[i] = mt.item[i];
		}
	}
	Matrix& operator=(const Matrix& mt)//重载赋值运算符
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
	Matrix operator*(const double x) const//函数重载，当*右侧是double时调用，矩阵各元素同乘该值
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
	friend Matrix operator*(const double x,const Matrix& mt)//友元函数，当矩阵*double时调用
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
	Matrix operator/(const double x) const //当矩阵/double时调用
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
	// set函数，用于接收二维数组作为参数进行赋值。利用 C++ 的模板特性，让编译器根据传入的实际二维数组类型来自动推导数组的维度信息
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

	// 重载set函数，用于接收初始化列表作为参数进行赋值
	 // 使用初始化列表赋值，这一段交给gpt了太麻烦了
	void set(std::initializer_list<std::initializer_list<double>> initList)
	{
		// 获取初始化列表的行和列
		row = initList.size();
		col = initList.begin()->size();

		// 确保矩阵大小和初始化列表匹配
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
		//返回矩阵的转置，不改变对象本身
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
		//交换对象矩阵的第x行和第y行，x和y从0开始计数
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
		//返回矩阵变换为的下三角形式，对象本身不变
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
						R.swap_row(j, j + k);//往下挨个找，直到确保(j,j)不为0
						change = true;
						break;
					}
				}
			}
			else
				change = true;
			if (!change)//这一列都找不到不为0的元素，直接处理下一列
				continue;
			for (int i = j + 1;i < row;i++)
			{
				//处理下方每一行
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
		//返回对象矩阵的伴随矩阵，对象本身不变
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
					//构造Temp矩阵，内容为this除去第i行第j列
					for (int m = 0;m < row;m++)
					{
						for (int n = 0;n < col;n++)
						{
							if (m == i || n == j)
								continue;
							//复制this的第i行第j列到Temp的第m>i?(m-1):m行第n>j?(n-1):n列
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
		//返回对象的行列式值，采用友元函数
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
						R.swap_row(j, j + k);//往下挨个找，直到确保(j,j)不为0
						change = true;
						row_change = !row_change;
						break;
					}
				}
			}
			else
				change = true;
			if (!change)//这一列都找不到不为0的元素，直接处理下一列
				continue;
			for (int i = j + 1;i < mt.row;i++)
			{
				//处理下方每一行
				double coe = R.item[i * mt.col + j] / R.item[j * mt.col + j];
				for (int k = 0;k < mt.col;k++)
				{
					R.item[i * mt.col + k] -= coe * R.item[j * mt.col + k];
				}
			}
		}
		//这里的R就是对象矩阵的下三角形式，与standard_foramt()函数不同的是这里需要保存采取交换两行操作的次数，每交换两行一次，行列式值取一次负。
		double results = 1.0;
		for (int i = 0;i < mt.row;i++)
		{
			results *= R.item[i * mt.col + i];//对角元素累乘得到行列式值
		}
		//在化下三角过程中，每进行一次交换两行的操作，行列式要取一次负号
		//运算符？：在这里的含义是：当row_change为真时返回results；当row_change为假时返回-results
		return row_change ? results : -results;
	}
	friend Matrix inv(const Matrix& mt)
	{
		//返回逆矩阵，采用友元函数形式
		if (mt.row != mt.col)
			throw std::range_error("Only square matrix has inverse matrix");
		//逆矩阵等于伴随矩阵除以行列式，运算符 / 已经重载，可以直接使用
		return mt.accompany() / det(mt);
	}

};

#endif