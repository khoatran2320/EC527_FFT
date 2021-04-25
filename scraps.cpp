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