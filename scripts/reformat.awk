BEGIN{
   print_on = 1;
}

($1=="/*__BEGIN_EXPERIMENTAL_SECTION__*/"){
   print_on = 0;
}

($1=="#__BEGIN_EXPERIMENTAL_SECTION__#"){
   print_on = 0;
}

($1=="/*___END_EXPERIMENTAL_SECTION___*/"){
   print_on = 1;
}

($1=="#___END_EXPERIMENTAL_SECTION___#"){
   print_on = 1;
}

($1=="/*UNCOMMENT"){
   getline;
   getline;
   while ($1 != "#endif"){
      print;
      getline;
   }
   getline;
}

($1=="/*" && $8=="Common"){
  getline;
  printf("/* This software is licensed under the Common Public License Version 1.0.    */\n");
}

($1=="SYMPHONYROOT" && $3=="${HOME}/SYMPHONY"){
   getline;
   printf("SYMPHONYROOT = ${HOME}/SYMPHONY-5.0\n");
}

($1=="USERROOT" && $3=="${SYMPHONYROOT}/Template"){
   getline;
   printf("USERROOT = ${SYMPHONYROOT}/USER\n");
}
     
($1=="USERROOT" && $3=="${SYMPHONYROOT}/Vrp"){
   getline;
   printf("USERROOT = ${SYMPHONYROOT}/VRP-5.0\n");
}

($1=="USERROOT" && $3=="${SYMPHONYROOT}/CNRP"){
   getline;
   printf("USERROOT = ${SYMPHONYROOT}/CNRP-5.0\n");
}

($1=="USERROOT" && $3=="${SYMPHONYROOT}/MCKP"){
   getline;
   printf("USERROOT = ${SYMPHONYROOT}/MCKP-5.0\n");
}
     
($1 == "USERROOT" && $3=="${SYMPHONYROOT}/SPP"){
   getline;
   printf("USERROOT = ${SYMPHONYROOT}/SPP-5.0\n");
}
     
($1=="USERROOT" && $3=="${SYMPHONYROOT}/SPP+CUTS"){
   getline;
   printf("USERROOT = ${SYMPHONYROOT}/SPP+CUTS-5.0\n");
}
     
($1=="USERROOT" && $3=="${SYMPHONYROOT}/MATCH"){
   getline;
   printf("USERROOT = ${SYMPHONYROOT}/MATCH-5.0\n");
}
     
($1=="USERROOT" && $3=="${SYMPHONYROOT}/MPP"){
   getline;
   printf("USERROOT = ${SYMPHONYROOT}/MPP-5.0\n");
}
     
($1=="/*" && $2=="accompanying"){
  getline;
  printf("/* Please see accompanying file for terms.                                   */\n");
}

($1!="/*___END_EXPERIMENTAL_SECTION___*/" && $1!="/*__BEGIN_EXPERIMENTAL_SECTION__*/" && $1!="/*UNCOMMENT*/" && $1!="#___END_EXPERIMENTAL_SECTION___#" &&
$1!="#___END_EXPERIMENTAL_SECTION___#"){
   if (print_on){
      print;
   }
}

END{
}







