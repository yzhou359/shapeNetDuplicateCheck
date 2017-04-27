#ifndef MATRIX2D_H
#define MATRIX2D_H

template<class T> class Matrix2D {

private:

	T* m;
	int rows, cols;

public:

	Matrix2D() {
		m = 0;
		rows = 0;
		cols = 0;
	}

	Matrix2D(T* m, int rows, int cols) {
		this->m = m;
		this->rows = rows;
		this->cols = cols;
	}

	inline T operator()(int i, int j) const{
        
        if (i < 0 || i >= rows || j < 0 || j >= cols)
            mexErrMsgTxt("Matrix2D index out of bounds.");
        
		return m[i + j*rows];
	}

	inline T& operator()(int i, int j) {
        
        if (i < 0 || i >= rows || j < 0 || j >= cols)
            mexErrMsgTxt("Matrix2D index out of bounds.");
        
		return m[i + j*rows];
	}
	
	inline int getRows() const {
		return rows;
	}

	inline int getCols() const {
		return cols;
	}

};

#endif