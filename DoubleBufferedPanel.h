#pragma once

public ref class DoubleBufferedPanel : System::Windows::Forms::Panel{
   public:
      void SetStyle(System::Windows::Forms::ControlStyles flag, System::Boolean value){
         System::Windows::Forms::Control::SetStyle(flag, value);
      }
};

