.PHONY: clean All

All:
	@echo "----------Building project:[ HW3NoGUI_Integrals - Debug ]----------"
	@cd "HW3NoGUI_Integrals" && $(MAKE) -f  "HW3NoGUI_Integrals.mk"
clean:
	@echo "----------Cleaning project:[ HW3NoGUI_Integrals - Debug ]----------"
	@cd "HW3NoGUI_Integrals" && $(MAKE) -f  "HW3NoGUI_Integrals.mk" clean
