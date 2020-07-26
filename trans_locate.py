
import constant as const
time =   0   #sdf_locate
locate = (time*const.dt_snapshot - const.window_start_time)*3e8*1e6

print(locate)
locate = 0 #locate
t = (locate/3e8 + const.window_start_time)/const.dt_snapshot
print(t)
