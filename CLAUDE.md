# CLAUDE.md

## Project Overview

FELsim is a Free Electron Laser simulation application with a Python backend (FastAPI/Uvicorn) and a frontend (Vite).

## Project Structure

- `backend/` - Python backend with FastAPI
- `fel-app/` - Frontend application (Vite)
- `fields/` - Field data files
- `beam_excel/` - Beam configuration Excel files
- `manuals/` - Reference documentation

## Running Locally

Backend:
```bash
cd backend/
pip install -r requirements.txt
uvicorn felAPI:app --host=127.0.0.1 --port=8000 --reload
```

Frontend:
```bash
cd fel-app/
npm install
npm run dev
```

## Reference Documentation

- [RF-Track Reference Manual](manuals/RF_Track_reference_manual.pdf) - Documentation for RF-Track beam dynamics simulation

## CTest-Like Orchestration

For sustained development tasks, define CTest-style acceptance criteria before coding. Work in a loop: run tests, fix failures, re-test until green. Applies to FELsim as a whole and to subcomponents (backend API, frontend, transport line simulation, field handling, etc.).

Example acceptance tests for this project:
- **Deterministic**: backend starts (`uvicorn felAPI:app`), frontend builds (`npm run build`), `pip install -r requirements.txt` succeeds
- **Output validation**: simulation output matches reference beam parameters within tolerance
- **AI-judged**: Haiku evaluates API endpoint correctness or physics model implementation (via `ai_judge.py` in `~/.claude-orchestrator/`)

See global CLAUDE.md for the full pattern description.

## Coding Style Guidelines

- Write code aligned with the standards of a senior software engineer or senior computational physicist: favour clarity, maintainability, and best practices over cleverness
- Avoid patterns characteristic of AI-generated code (e.g., excessive comments explaining obvious operations, overly verbose variable names, unnecessary abstractions, repetitive boilerplate)
- Use "Eremey Valetov" for authorship attribution in file headers, not placeholders like "[Your Name]"
- Follow existing project conventions and style when extending the codebase
- Prefer concise, self-documenting code with comments reserved for non-obvious logic
