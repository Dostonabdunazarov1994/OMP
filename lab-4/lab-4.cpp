#include <iostream>
#include <string>
#include <math.h>
#include <time.h>
#include <omp.h>
using namespace std;

void SingleTask1(double** A, int N);
void ParallelforTask1(double** A, int N, int k);
void ReductionTask1(double** A, int N, int k);
void AtomicTask1(double** A, int N, int k);
void CriticalTask1(double** A, int N, int k);
double** ArrGenerator(const int N);
void array_destroyer(double** A, const int N);
void MultiplyMatrix(double** A, double** B, double** C, int N, int k);
/*
Задание 1.
Работа происходит с квадратной матрицей А вещественных чисел размером NxN.Матрица вычисляется по формуле.
A_ij = Sin(i) + Cos(j), i, j = 0, …, N - 1
Найти среднее арифметическое элементов матрицы.
Порядок выполнения :
1. В качестве N выбрать максимально доступное для реализации число, кратное 100.
2. Реализовать последовательный вариант и определить время его работы.
3. Выделить участки последовательного и распараллеливаемого кода.
4. Выполнить распараллеливание в OpenMP(четырьмя способами : parallel for, reduction, atomic, критические области)
   и оценить их эффективность(на основе среднего из 20 экспериментов).Количество потоков М = 2, 4, 8.
5. В качестве контроля корректности распараллеливания сравнить найденные средние с полученным в последовательном
   режиме.
*/
//Задание 2.
//Проведите серию экспериментов на персональном компьютере по исследованию масшта - бируемости 
//OpenMP - программ умножения квадратных матриц размера NхN

int main()
{
	setlocale(LC_ALL, "RUS");
	const int N = 5000;
	int M[3] = { 2, 4, 8 };
	double** A = ArrGenerator(N); //матрица NxN
	cout << "single:  " << endl;
	SingleTask1(A, N); //последовательное выполнения
	cout << "------------------------------------------" << endl;
	cout << endl << "parallel for: " << endl;
	for (int i = 0; i < 3; i++) {//параллельное выполнения
		cout << "при k = " << M[i] << " ";
		ParallelforTask1(A, N, M[i]);
	}
	cout << "------------------------------------------" << endl;
	cout << endl << "reduction: " << endl; // reduction
	for (int i = 0; i < 3; i++) {
		cout << "при k = " << M[i] << " ";
		ReductionTask1(A, N, M[i]);
	}
	cout << "------------------------------------------" << endl;
	cout << endl << "atomic: " << endl;
	for (int i = 0; i < 3; i++) { //atomic
		cout << "при k = " << M[i] << " ";
		AtomicTask1(A, N, M[i]);
	}
	cout << "------------------------------------------" << endl;
	cout << endl << "critical: " << endl;
	for (int i = 0; i < 3; i++) {//critical
		cout << "при k = " << M[i] << " ";
		CriticalTask1(A, N, M[i]); 
	}
	cout << "------------------------------------------" << endl;
	array_destroyer(A, N); //деинициализация

	int N1[4] = { 100, 1000, 5000, 7000 }; //кол-во элементов матриц 
	int k[9] = { 1, 2, 4, 8, 10, 16, 20, 24, 30 }; //кол-во потоков
	double** A1, ** B, ** C; //матрицы для умножения 
	int n = 0;
	for (int i = 0; i < (sizeof(N1) / sizeof(N1[0])); i++) {
		A1 = ArrGenerator(N1[i]); //инициализация
		B = ArrGenerator(N1[i]);
		C = ArrGenerator(N1[i]);
		n = N1[i];
		cout << "N = " << N1[i] << endl;
		for (int j = 0; j < (sizeof(k) / sizeof(k[0])); j++) {
			cout << "   k = " << k[j] << " ";
			MultiplyMatrix(A1, B, C, N1[i], k[j]); //функция для умножения матриц
		}
		cout << endl;
	}	
	cout << "end";
	array_destroyer(A1, n); //деинициализация
	array_destroyer(B, n);
	array_destroyer(C, n);

	return 0;
}

double** ArrGenerator(const int N) { //функция для иницилизации матриц
	double** arr = new double* [N];
	for (int i = 0; i < N; i++) {
		arr[i] = new double[N];
	}
	return arr;
}
void array_destroyer(double** arr, const int N) { //функция для удаления матриц
	for (int i = 0; i < N; i++) {
		delete[] arr[i];
	}
	delete[] arr;
}
void SingleTask1(double** A, int N) { //последовательное выполнения
	double sum = 0.0;
	int count = 0;
	clock_t start = clock();
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			A[i][j] = sin(i) + cos(j);
			sum += A[i][j];
			count++;
		}
	}
	clock_t end = clock();
	double dt = (double)(end - start) / CLOCKS_PER_SEC;
	cout << "время: " << dt << "cек;  Ср. арифметическое: " << sum / count << endl;
}

void ParallelforTask1(double** A, int N, int k) {//параллел c помошью parallel for
	double sum = 0.0;
	clock_t start = clock();
#pragma omp parallel for num_threads(k)
	for (int i = 0; i < N; i++) {
		double sumStr = 0.0;
		for (int j = 0; j < N; j++) {
			A[i][j] = sin(i) + cos(j);
			sumStr += A[i][j];
		}
		sum += sumStr;
	}
	clock_t end = clock();
	double dt = (double)(end - start) / CLOCKS_PER_SEC;
	cout << "время: " << dt << "cек;  Ср. арифметическое: " << sum / (N * N) << endl;
}

void ReductionTask1(double** A, int N, int k) {//c помошью reduction
	double sum = 0.0;
	clock_t start = clock();
#pragma omp parallel for num_threads(k) reduction(+ : sum)
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			A[i][j] = sin(i) + cos(j);
			sum += A[i][j];
		}
	}
	clock_t end = clock();
	double dt = (double)(end - start) / CLOCKS_PER_SEC;
	cout << "время: " << dt << "cек;  Ср. арифметическое: " << sum / (N * N) << endl;
}

void AtomicTask1(double** A, int N, int k) {//c помошью atomic
	double sum = 0;
	clock_t start = clock();
#pragma omp parallel for num_threads(k)
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			A[i][j] = sin(i) + cos(j);
#pragma omp atomic
			sum += A[i][j];
		}
	}
	clock_t end = clock();
	double dt = (double)(end - start) / CLOCKS_PER_SEC;
	cout << "время: " << dt << "cек;  Ср. арифметическое: " << sum / (N * N) << endl;
}

void CriticalTask1(double** A, int N, int k) {//c помошью critical
	double sum = 0.0;
	clock_t start = clock();
#pragma omp parallel for num_threads(k)
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			A[i][j] = sin(i) + cos(j);
#pragma omp critical
			sum += A[i][j];
		}
	}
	clock_t end = clock();
	double dt = (double)(end - start) / CLOCKS_PER_SEC;
	cout << "время: " << dt << "cек;  Ср. арифметическое: " << sum / (N * N) << endl;
}

void MultiplyMatrix(double** A, double** B, double** C, int N, int k) { //умножения матриц
	clock_t start = clock();
#pragma omp parallel for num_threads(k) schedule(guided)
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			A[i][j] = rand() % 11;
			B[i][j] = rand() % 11;
		}
	}
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			C[i][j] = A[i][j] * B[j][i];
		}
	}
	clock_t end = clock();
	double dt = (double)(end - start) / CLOCKS_PER_SEC;
	cout << "время: " << dt << " сек" << endl;
}