#include "SimWindow.h"

using namespace System;
using namespace System::Windows::Forms;
#include<iostream>
#include<Eigen\Sparse>
#include<math.h>
#include<typeinfo>
#include<ctime>
#include "myutilities.h"
#include"DoubleBufferedPanel.h"
double f(array<double, 2>^);
[STAThreadAttribute]
void Main(array<String^>^ args) {
	Application::EnableVisualStyles();
	Application::SetCompatibleTextRenderingDefault(false);
	Sim::SimWindow SimWindow;
	Application::Run(%SimWindow);
   array<double, 2>^ a = gcnew array<double, 2>(10, 10);
   a[0, 0] = 5;
   std::cout << f(a);
   int x = 0;
}
double f(array<double, 2>^ b){
   return b[0, 0] + 1;
}