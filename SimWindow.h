#pragma once
#include<iostream>
#include<Eigen\Sparse>
#include<math.h>
#include<typeinfo>
#include<ctime>
#include "myutilities.h"
#include"DoubleBufferedPanel.h"
namespace Sim{
   
	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;
   using namespace System::Collections::Generic;
   #define SWAP(x,y) {array<double,2>^ tmp=x;x=y;y=tmp;}
	/// <summary>
	/// Summary for SimWindow
	/// </summary>
	public ref class SimWindow: public System::Windows::Forms::Form{
		public:
		   //SimWindow Constructor
		   SimWindow(void){
			   InitializeComponent();
            for (int i = 0; i < SIZE; i++){
               for (int j = 0; j < SIZE; j++){
                  _u[i,j] = 0.0;
                  _v[i,j] = 0.0;
                  _dens[i,j]=0.0;
                  _u_prev[i,j] = 0.0;
                  _v_prev[i,j] = 0.0;
                  _dens_prev[i,j] = 0.0;
               }
            }

            this->textBox1->Text += _diff;
            this->textBox2->Text += _dt;
            this->textBox3->Text += N;
            this->textBox4->Text += _viscosity;
			   
			   
			   _defaultimage = gcnew Bitmap(SIZE, SIZE);
            static_cast<DoubleBufferedPanel^>(SimPanel)->SetStyle(ControlStyles::AllPaintingInWmPaint, true);
            static_cast<DoubleBufferedPanel^>(SimPanel)->SetStyle(ControlStyles::DoubleBuffer, true);
            Color c = Color::FromArgb(0,119,190);
			   for (int Xcount = 0; Xcount < _defaultimage->Width; Xcount++){
				   for (int Ycount = 0; Ycount < _defaultimage->Height; Ycount++){
                  
					   _defaultimage->SetPixel(Xcount, Ycount, c);
				   }
			   }
			   SimPanel->BackgroundImage = _defaultimage;
			   timer->Start();
		   }
		protected:
		   /// <summary>
		   /// Clean up any resources being used.
		   /// </summary>
		   ~SimWindow(){
			   if (components){
				   delete components;
			   }
		   }
		private:
         System::Windows::Forms::Panel^  SimPanel;
      private: System::Windows::Forms::TextBox^  textBox1;
      private: System::Windows::Forms::TextBox^  textBox2;
      private: System::Windows::Forms::TextBox^  textBox3;
      private: System::Windows::Forms::TextBox^  textBox4;
               System::ComponentModel::IContainer^  components;
		   #pragma region Windows Form Designer generated code
		   /// <summary>
		   /// Required method for Designer support - do not modify
		   /// the contents of this method with the code editor.
		   /// </summary>
		   void InitializeComponent(void){
            this->components = (gcnew System::ComponentModel::Container());
            this->SimPanel = (gcnew System::Windows::Forms::Panel());
            this->timer = (gcnew System::Windows::Forms::Timer(this->components));
            this->textBox1 = (gcnew System::Windows::Forms::TextBox());
            this->textBox2 = (gcnew System::Windows::Forms::TextBox());
            this->textBox3 = (gcnew System::Windows::Forms::TextBox());
            this->textBox4 = (gcnew System::Windows::Forms::TextBox());
            this->SuspendLayout();
            // 
            // SimPanel
            // 
            this->SimPanel->BackgroundImageLayout = System::Windows::Forms::ImageLayout::Stretch;
            this->SimPanel->Location = System::Drawing::Point(10, 10);
            this->SimPanel->Name = L"SimPanel";
            this->SimPanel->Size = System::Drawing::Size(400, 400);
            this->SimPanel->TabIndex = 0;
            this->SimPanel->MouseClick += gcnew System::Windows::Forms::MouseEventHandler(this, &SimWindow::panel1_OnMouseClick);
            // 
            // timer
            // 
            this->timer->Interval = 10;
            this->timer->Tick += gcnew System::EventHandler(this, &SimWindow::timer_Tick);
            // 
            // textBox1
            // 
            this->textBox1->Location = System::Drawing::Point(416, 12);
            this->textBox1->Name = L"textBox1";
            this->textBox1->ReadOnly = true;
            this->textBox1->Size = System::Drawing::Size(106, 20);
            this->textBox1->TabIndex = 1;
            this->textBox1->Text = L"_diff = ";
            // 
            // textBox2
            // 
            this->textBox2->Location = System::Drawing::Point(416, 38);
            this->textBox2->Name = L"textBox2";
            this->textBox2->ReadOnly = true;
            this->textBox2->Size = System::Drawing::Size(106, 20);
            this->textBox2->TabIndex = 2;
            this->textBox2->Text = L"_dt = ";
            // 
            // textBox3
            // 
            this->textBox3->Location = System::Drawing::Point(416, 64);
            this->textBox3->Name = L"textBox3";
            this->textBox3->ReadOnly = true;
            this->textBox3->Size = System::Drawing::Size(106, 20);
            this->textBox3->TabIndex = 3;
            this->textBox3->Text = L"N = ";
            // 
            // textBox4
            // 
            this->textBox4->Location = System::Drawing::Point(416, 90);
            this->textBox4->Name = L"textBox4";
            this->textBox4->ReadOnly = true;
            this->textBox4->Size = System::Drawing::Size(106, 20);
            this->textBox4->TabIndex = 4;
            this->textBox4->Text = L"_visc = ";
            // 
            // SimWindow
            // 
            this->ClientSize = System::Drawing::Size(534, 421);
            this->Controls->Add(this->textBox4);
            this->Controls->Add(this->textBox3);
            this->Controls->Add(this->textBox2);
            this->Controls->Add(this->textBox1);
            this->Controls->Add(this->SimPanel);
            this->DoubleBuffered = true;
            this->Name = L"SimWindow";
            this->Text = L"SimWindow";
            this->ResumeLayout(false);
            this->PerformLayout();

         }
         ///
         ///Used in simulator.
         ///
         const int N = 200;
         const int SIZE = 202;
         const double _dt = .1;
         double _viscosity = 1;
         double _diff = 10000;
         array<double,2>^ _u         = gcnew array<double,2>(202,202);
         array<double,2>^ _u_prev    = gcnew array<double,2>(202,202);
         array<double,2>^ _v         = gcnew array<double,2>(202,202);
         array<double,2>^ _v_prev    = gcnew array<double,2>(202,202);
         array<double,2>^ _dens      = gcnew array<double,2>(202,202);
         array<double,2>^ _dens_prev = gcnew array<double,2>(202,202);
         System::Windows::Forms::Timer^  timer;
         Bitmap^ _defaultimage;
         int _dir = 1;

		#pragma endregion
		   void AddMouseEffect(Point location){
			   //Color c = GetRandomColor();
			   //for (int Xcount = 0; Xcount < _defaultimage->Width; Xcount++){
			   //	for (int Ycount = 0; Ycount < _defaultimage->Height; Ycount++){
			   //		this->_defaultimage->SetPixel(Xcount, Ycount, c);
			   //	}
			   //}
			   int i = (int)(Math::Round(location.X*((double)N/(double)this->SimPanel->Width)) + 1);
			   int j = (int)(Math::Round(location.Y*((double)N/(double)this->SimPanel->Height)) + 1);
            for(int k = Math::Max(0,i - 5); k<= i+5; k++){
               for(int l = Math::Max(0,j-5); l<= j+5; l++){
                  _dens_prev[k,l] += _diff;
                  _v_prev[k,l] +=_diff*_dir;
                  _u_prev[k,l]+=_diff*_dir;
                }
            }
			   _dir*= -1;
		   }
		   System::Void panel1_OnMouseClick(System::Object^ sender, System::Windows::Forms::MouseEventArgs^ e){
			   AddMouseEffect(e->Location);
		   }
		   System::Void timer_Tick(System::Object^  sender, System::EventArgs^  e){
               Bitmap^ bmp = gcnew Bitmap(_defaultimage);
            vel_step(N,_u,_v,_u_prev,_v_prev,_viscosity,_dt);
            dens_step(N,_dens,_dens_prev,_u,_v,_diff,_dt);

            for (int i = 1; i <= N; i++){
               for (int j = 1; j <= N; j++){
                  double density = _dens[i,j];
                  double velocity = _v[i,j];
                  if (abs(density) >= 0.5 || abs(velocity) >= 0.5){
                     int r = Math::Min(255, Math::Max(0, (int)Math::Round(Math::Abs(10*velocity))));
                     int g = Math::Min(255, Math::Max(0, 119 + (int)Math::Round(Math::Abs(.5*velocity)) + (int)Math::Round(Math::Abs(.5*density))));
                     int b = Math::Min(255, Math::Max(0, 190 + (int)Math::Round(density)));
                     
                     bmp->SetPixel(i - 1, j - 1, Color::FromArgb(r, g, b));
                  }

               }
            }
            Bitmap^ bmpt = bmp;
            SimPanel->BackgroundImage = bmpt;
            //if(bmp != nullptr){ delete bmp; }
         }
         void add_source(int N, cli::array<double, 2>^ x, cli::array<double, 2>^ s, double dt){
            int i, size = 49;
            for (i = 0; i<N + 2; i++){
               for (int j = 0; j < N + 2; j++){
                  x[i, j] += dt * s[i, j];
               }

            }
         }

         void diffuse(int N, int b, cli::array<double, 2>^ x, cli::array<double, 2>^ x0, double diff, double dt){
            int i, j, k;
            double a = dt * diff*N*N;
            for (k = 0; k<20; k++){
               for (i = 1; i <= N; i++){
                  for (j = 1; j <= N; j++){
                     x[i, j] = (x0[i,j] + a * (x[i - 1, j] + x[i + 1, j] + x[i, j - 1] + x[i, j + 1])) / (1 + 4 * a);

                  }
               }
            }
            set_bnd(N, b, x);
         }

         void advect(int N, int b, cli::array<double, 2>^ d, cli::array<double, 2>^ d0, cli::array<double, 2>^ u, cli::array<double, 2>^ v, double dt){
            int i, j, i0, j0, i1, j1;
            double x, y, s0, t0, s1, t1, dt0;
            dt0 = dt * N;
            for (i = 1; i <= N; i++){
               for (j = 1; j <= N; j++){
                  x = i - dt0 * u[i, j];
                  y = j - dt0 * v[i, j];
                  if (x<.5) x = .5;
                  if (x>N + .5) x = N + .5;
                  i0 = (int)x;
                  i1 = i0 + 1;
                  if (y < 0.5) y = 0.5;
                  if (y > N + 0.5) y = N + 0.5;
                  j0 = (int)y;
                  j1 = j0 + 1;
                  s1 = x - i0;
                  s0 = 1 - s1;
                  t1 = y - j0;
                  t0 = 1 - t1;
                  d[i, j] = s0 * (t0*d0[i0, j0] + t1 * d0[i0, j1])
                     + s1 * (t0*d0[i1, j0] + t1 * d0[i0, j1]);

               }
            }
            set_bnd(N, b, d);
         }

         //x: dens, x0: dens_prev, u: u, v: v
         void dens_step(int N, cli::array<double, 2>^ x, cli::array<double, 2>^ x0, cli::array<double, 2>^ u, cli::array<double, 2>^ v, double diff, double dt){
            add_source(N, x, x0, dt);
            SWAP(x0, x);
            diffuse(N, 0, x, x0, diff, dt);
            SWAP(x0, x);
            advect(N, 0, x, x0, u, v, dt);

         }
         //u, v, u_prev, v_prev
         void vel_step(int N, cli::array<double, 2>^ u, cli::array<double, 2>^ v, cli::array<double, 2>^ u0, cli::array<double, 2>^ v0, float visc, double dt){
            add_source(N, u, u0, dt);
            add_source(N, v, v0, dt);
            SWAP(u0, u);
            diffuse(N, 1, u, u0, visc, dt);
            SWAP(v0, v);
            diffuse(N, 2, v, v0, visc, dt);
            project(N, u, v, u0, v0);
            SWAP(u0, u);
            SWAP(v0, v);
            advect(N, 1, u, u0, u0, v0, dt);
            advect(N, 2, v, v0, u0, v0, dt);
            project(N,u,v,u0,v0);
         }

         void project(int N, cli::array<double, 2>^ u, cli::array<double, 2>^ v, cli::array<double, 2>^ p, cli::array<double, 2>^ div){
            int i, j, k;
            double h;
            h = 1.0 / N;
            for (i = 1; i <= N; i++){
               for (j = 1; j <= N; j++){
                  div[i, j] = -0.5*h*(u[i + 1, j] - u[i - 1, j]
                                      + v[i, j + 1] - v[i, j - 1]);
                  p[i, j] = 0;
               }
            }
            set_bnd(N, 0, div);
            set_bnd(N, 0, p);
            for (k = 0; k<20; k++){
               for (i = 1; i <= N; i++){
                  for (j = 1; j <= N; j++){
                     p[i, j] = (div[i, j] + p[i - 1, j] + p[i + 1, j]
                                + p[i, j - 1] + p[i, j + 1]) / 4;
                  }
               }
               set_bnd(N, 0, p);
            }
            for (i = 1; i <= N; i++){
               for (j = 1; j <= N; j++){
                  u[i, j] -= 0.5*(p[i + 1, j] - p[i - 1, j]) / h;
                  u[i, j] -= 0.5*(p[i, j + 1] - p[i, j - 1]) / h;
               }
            }
            set_bnd(N, 1, u);
            set_bnd(N, 2, v);
         }

         void set_bnd(int N, int b, cli::array<double, 2>^ x){
            for (int i = 1; i <= N; i++){
               x[0, i] = b == 1 ? -x[1, i] : x[1, i];
               x[N + 1, i] = b == 1 ? -x[N, i] : x[N, i];
               x[i, 0] = b == 2 ? -x[i, 1] : x[i, 1];
               x[i, N + 1] = b == 2 ? -x[i, N] : x[i, N];
            }
            x[0, 0] = 0.5 * (x[1, 0] + x[0, 1]);
            x[0, N + 1] = 0.5 * (x[1, N + 1] + x[0, N]);
            x[N + 1, 0] = 0.5 * (x[N, 0] + x[N + 1, 1]);
            x[N + 1, N + 1] = 0.5 * (x[N, N + 1] + x[N + 1, N]);
         }
	};
}

