import platform

os_info = open('data/operating-system.txt', 'w')
os_name = platform.system()

if os_name == "Linux":
    os_info.write("Linux")
    
elif os_name == "Darwin":
     os_info.write("macOS")
     
elif os_name == "Windows":
     os_info.write("Windows")
     
else:
    print("Unknown OS")

os_info.close()