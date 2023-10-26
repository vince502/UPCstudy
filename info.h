bool pi(float x, float y){
    if(y>0){
        if( (log(y) + 3*x < 2.2 ) && (log(y) + 7.5*x > 1.4) ) return true;
    }
	if(y<2.8){
		if( x > 0.4) return true;
	}
    return false;
};
bool ka(float x, float y){
    if(y>3){
        if( (log(y) + 3*x > 2.5 ) && (log(y) + 2.64*x < 2.9) ) return true;
    }
    return false;
};
bool pr(float x, float y){
    if(y>2.65){
        if( (log(y) + 2.64*x > 3.03) ) return true;
    }
    return false;
};
bool e(float x, float y){
	return !(pi(x,y) || ka(x,y) || pr(x,y));
};
