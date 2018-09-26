import sys

def print(max_val, max_bar_val, cur_val):
#    print.old_num_asts += 1
    bar = ""
    wheel = {0: "/", 1: "-", 2: "\\", 3: "|", 4: "/", 5: "-", 6: "\\", 7: "|"}
    num_ast = int(max_bar_val * cur_val / max_val)
    for j in range(0, num_ast-1):
        bar += "*"
    bar += wheel[cur_val%8]
    for j in range(num_ast, max_bar_val + 1):
        bar += "-"
        #        print(bar, end='\r', flush=True)
    sys.stdout.write("\r" + bar)
    sys.stdout.flush()

#print.old_num_asts = 0