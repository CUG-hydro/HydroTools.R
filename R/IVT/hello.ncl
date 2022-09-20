
load "main_ivt.ncl"

hello_world()

x = True

list_vars

if isvar("x") then
  print(x)
end if 

param = /
; x = 2
; varname = "x"
; y="$x$"
; printVarSummary(y)
; print(y)
; print(y)

; vapor_pres_rh(rh, esat)


; undef("hello")
; function hello(x)

; begin
; if ismissing(x) then
;   x = 2
; end if 

; print(x)
; ; local local_variables ; optional, but recommended
; ; statements
; return(x)
; end

; """
; This document
; """

; if .not. isdefined("x") then
;   print(1)
; else
;   print(x)
; end if

; ; print(x)
; ; print(y)
; ; y = hello(missing)

; ; undef("function_name")
; ; function function_name(declaration_list)
; ; local local_variables ; optional, but recommended
; ; begin
; ; statements
; ; return(return_variable)
; ; end
