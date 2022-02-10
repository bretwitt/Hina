#ifndef _FUNCTION_H
#define _FUNCTION_H

typedef double (*function_ptr)(double,void*);

class Function {

private:
    function_ptr funct;

public:

    Function(function_ptr funct) {
        this->funct = funct;
    };

    Function() {
        this->funct = nullptr;
    };

    virtual function_ptr getFunction() {
        return funct;
    };

};

#endif