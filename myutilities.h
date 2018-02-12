#pragma once
#define COLORBASE (256)
static System::Drawing::Color GetRandomColor(void){
	srand(time(0));
	int red, green, blue;
	red = rand() % COLORBASE;
	green = rand() % COLORBASE;
	blue = rand() % COLORBASE;

	return System::Drawing::Color::FromArgb(red, green, blue);
}
