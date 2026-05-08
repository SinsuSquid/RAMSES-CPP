import os
import sys
import argparse
import textwrap

def port_with_gemini(fortran_code, system_prompt):
    try:
        import google.generativeai as genai
    except ImportError:
        print("Error: google-generativeai package not found. Run 'pip install google-generativeai'")
        sys.exit(1)

    api_key = os.getenv("GEMINI_API_KEY")
    if not api_key:
        print("Error: GEMINI_API_KEY environment variable not set.")
        sys.exit(1)

    genai.configure(api_key=api_key)
    model = genai.GenerativeModel('gemini-1.5-pro')
    
    full_prompt = f"{system_prompt}\n\nTranslate the following Fortran code into the C++ override class format specified above:\n\n{fortran_code}"
    
    response = model.generate_content(full_prompt)
    return response.text

def port_with_openai(fortran_code, system_prompt):
    try:
        from openai import OpenAI
    except ImportError:
        print("Error: openai package not found. Run 'pip install openai'")
        sys.exit(1)

    api_key = os.getenv("OPENAI_API_KEY")
    if not api_key:
        print("Error: OPENAI_API_KEY environment variable not set.")
        sys.exit(1)

    client = OpenAI(api_key=api_key)
    response = client.chat.completions.create(
        model="gpt-4-turbo",
        messages=[
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": f"Translate this Fortran code:\n\n{fortran_code}"}
        ]
    )
    return response.choices[0].message.content

def main():
    parser = argparse.ArgumentParser(description="RAMSES-CPP AI Patch Porter Tool ✨")
    parser.add_argument("input", help="Path to the legacy Fortran (.f90) patch file")
    parser.add_argument("--solver", required=True, help="Target solver name (e.g., Cooling, Hydro, Poisson)")
    parser.add_argument("--provider", default="gemini", choices=["gemini", "openai"], help="LLM provider (default: gemini)")
    
    args = parser.parse_args()

    # 1. Load Rosetta Stone Prompt
    script_dir = os.path.dirname(os.path.realpath(__file__))
    prompt_path = os.path.join(script_dir, "rosetta_prompt.txt")
    with open(prompt_path, "r") as f:
        system_prompt = f.read()

    # 2. Load Fortran Code
    with open(args.input, "r") as f:
        fortran_code = f.read()

    print(f"🚀 Porting {args.input} using {args.provider}...")

    # 3. Call LLM
    if args.provider == "gemini":
        cpp_code = port_with_gemini(fortran_code, system_prompt)
    else:
        cpp_code = port_with_openai(fortran_code, system_prompt)

    # 4. Clean up output (remove markdown code blocks)
    cpp_code = cpp_code.replace("```cpp", "").replace("```", "").strip()

    # 5. Save to cpp_patches/
    repo_root = os.path.abspath(os.path.join(script_dir, "..", ".."))
    out_dir = os.path.join(repo_root, "cpp_patches")
    os.makedirs(out_dir, exist_ok=True)
    
    out_file = os.path.join(out_dir, f"Patch{args.solver}Solver.cpp")
    with open(out_file, "w") as f:
        f.write(cpp_code)

    print(f"✨ Successfully generated patch: {out_file}")
    print("👉 Re-run CMake to compile the new patch!")

if __name__ == "__main__":
    main()
