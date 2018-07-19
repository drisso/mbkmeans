#ifndef __Week4_functions__
#define __Weel4_functions__

int type_integer_index(SEXP data);
int type_numeric_index(SEXP data);
std::string make_to_string(const Rcpp::RObject& str);
std::string get_class(const Rcpp::RObject& incoming);
Rcpp::RObject get_safe_slot(const Rcpp::RObject& incoming, const std::string& slotname);

#endif