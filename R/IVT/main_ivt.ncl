; """
; 如果`name` attribute 不存在，则设置为默认值`default`
; 
; ## Examples
; ```ncl
; p = True
; p@x = 1
; ; param@y = 2
; default_att(p, "y", 3)
; default_atts(p, [/"z"/], [/4/])
; print(p)
; ```
; """
undef("default_att")
procedure default_att(x, name:string, value)
begin
if(.not.isatt(x, name)) then
  x@$name$ = value
end if
end ; end default_atts


undef("default_atts")
procedure default_atts(x, names, values)
local name, value
begin
n = ListCount(names)

do i = 0, n-1
  name = names[i]
  value = values[i]
  ; print(i)
  ; print(name)
  if(.not.isatt(x, name)) then
    x@$name$ = value
  end if  
end do
end ; end default_atts


undef("print_var")
procedure print_var(x, name)
begin
  print("[" + name + "] " + "-----")
  printVarSummary(x)                          ; (time,level,lat,lon); m/s
  printMinMax(x,0)
end

y = 1
print_var(y, "y")
