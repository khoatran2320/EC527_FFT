int main(){
	int i,j;
	float * a = new (std::align_val_t(32)) float[3];//{1.1,2.2,3.3};
	for (i = 0; i < 3; i++){
		a[i * sizeof(float)] = 1.1;
		cout << a[i * sizeof(float)] << " ";
	}
	cout << endl;
	int freq = 100;
	int periods = 1;
	int len = 3;//int) (2 * PI * periods * freq) + 1;
	__m256 sin_x[len];
//	//for (i = 0; i < len; i++){
	//	float  f = 1.1;
		float* fp_ptr = a;
		_mm256_store_ps(a, sin_x[0]);
//	//}
	return 0;
}


int main(){
	int i,j;
	float * a = new (std::align_val_t(32)) float[3];//{1.1,2.2,3.3};
	float * out = new (std::align_val_t(32)) float[3];
	for (i = 0; i < 3; i++){
		a[i * sizeof(float)] = 3.14159/2;
		cout << a[i * sizeof(float)] << " ";
	}
	cout << endl;
	int freq = 100;
	int periods = 1;
	int len = 3;//int) (2 * PI * periods * freq) + 1;
	//__m256 sin_x[len];
//	//for (i = 0; i < len; i++){
	//	float  f = 1.1;
		const float * fp_ptr = a;
		__m256 sin_x = _mm256_load_ps(fp_ptr);
		//__m256 sin_y;
		__m256d sin_y = _ZGVcN4v_sin(sin_x);
		
//		_mm256_store_ps(out, sin_x);
	//	_mm256_load_ps(sin_x[0]);
		for (i = 0; i < 3; i++) {
		cout << out[0 + i * sizeof(float)] << " ";
	}
	cout << endl;

//	//}
	return 0;
}

void proof_of_concept(){
		int freq = 10;
	int periods = 1;
	int len = (int) 2 * PI * periods * freq;
	int index = 0;

	__m256 blocks[len];

	for (int j = 0; j < len; j += 8) {
	//	cout << "j : " << j << endl;

	float * a_float = new (std::align_val_t(32)) float[8];
	float * b_float = new (std::align_val_t(32)) float[8];

	for (int i = 0; i < 8; i++) a_float[i] = 10;// * sin((float)i/freq);
	for (int i = 0; i < 8; i++) b_float[i] = sin(((float)i + index)/freq);


//	for (int i = 0; i < len; i++) cout << a_float[i] << endl;


	__m256 a = _mm256_load_ps(a_float);
	__m256 b = _mm256_load_ps(b_float);
	__m256 c = _mm256_sub_ps(b,a);

	blocks[0] = c;
	float * c_float = new (std::align_val_t(32)) float[8];	


	_mm256_store_ps(c_float,c);

	//for (int i = 0; i < 8; i++) cout << c_float[i] << endl;
	index += 8;
}
	float * d_float = new (std::align_val_t(32)) float[8];	

	_mm256_store_ps(d_float, blocks[0]);
	for (int i = 0; i < 8; i++) cout << d_float[i] << endl;
}