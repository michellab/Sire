

def fixMB(mb):
    mb.add_declaration_code( "#include \"_IO_load.h" )
    mb.add_registration_code( "register_SireIO_load_function();")
